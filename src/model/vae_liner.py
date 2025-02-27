import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from torch.utils.data import DataLoader, TensorDataset
from torch.nn.utils.rnn import pad_sequence
import torch.nn.functional as F
import os
import torch.distributed as dist
from tqdm import tqdm


class VAEModel2(nn.Module):
    def __init__(self, sequence_len, input_dim, latent_dim=200, lstm_units=128):
        super(VAEModel2, self).__init__()
        self.sequence_len = sequence_len
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.lstm_units = lstm_units

        # ----------------- 编码器（双向LSTM） -----------------
        self.encoder_lstm = nn.LSTM(
            input_size=input_dim,
            hidden_size=lstm_units,
            batch_first=True,
            bidirectional=True  # 编码器双向
        )
        # 全连接层输入维度：双向拼接后的维度
        self.fc_mu = nn.Linear(2 * lstm_units, latent_dim)  # 2*lstm_units
        self.fc_log_var = nn.Linear(2 * lstm_units, latent_dim)

        # ----------------- 解码器（双向LSTM） -----------------
        self.decoder_lstm = nn.LSTM(
            input_size=latent_dim,
            hidden_size=lstm_units,
            batch_first=True,
            bidirectional=True  # 解码器双向（关键修改）
        )
        # 解码器全连接层输入维度：双向拼接后的维度
        self.decoder_fc = nn.Linear(2 * lstm_units, input_dim)  # 2*lstm_units

    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x, mask):
        # ----------------- 编码器部分 -----------------
        # 输入 x 形状: (batch_size, sequence_len, input_dim)
        lstm_out, _ = self.encoder_lstm(x)  # 输出形状: (batch_size, seq_len, 2*lstm_units)
        mu = self.fc_mu(lstm_out[:, -1, :])  # 取最后一个时间步的双向特征
        log_var = self.fc_log_var(lstm_out[:, -1, :])
        z = self.reparameterize(mu, log_var)  # 形状: (batch_size, latent_dim)

        # ----------------- 解码器部分 -----------------
        # 扩展潜在变量 z 为序列输入
        z = z.unsqueeze(1).repeat(1, self.sequence_len, 1)  # (batch_size, seq_len, latent_dim)
        lstm_out, _ = self.decoder_lstm(z)  # 输出形状: (batch_size, seq_len, 2*lstm_units)
        recon_x = self.decoder_fc(lstm_out)  # 形状: (batch_size, seq_len, input_dim)

        # 应用掩码过滤填充部分
        recon_x = recon_x * mask.unsqueeze(-1)
        return recon_x, mu, log_var

    def loss_function(self, recon_x, x, mu, log_var, mask):
        # 重构损失（仅计算有效位置）
        reconstruction_loss = F.mse_loss(recon_x, x, reduction='none')
        reconstruction_loss = torch.sum(reconstruction_loss * mask.unsqueeze(-1))
        # KL散度损失
        kl_divergence = -0.5 * torch.mean(1 + log_var - mu.pow(2) - log_var.exp())
        return reconstruction_loss + kl_divergence

    def generate_latent_features(self, x):
        # 生成双向LSTM的完整序列特征
        lstm_out, _ = self.encoder_lstm(x)  # 形状: (batch_size, seq_len, 2*lstm_units)
        mu = self.fc_mu(lstm_out[:, -1, :])  # 取最后一个时间步的双向特征
        return mu


class VAETrainer:
    def __init__(self, model, data, device,batch_size=128, epochs=50, learning_rate=1e-3,multiple_gpu=True):
        # init_dist()
        self.multiple_gpu=multiple_gpu
        if device.type != 'cpu' and torch.cuda.device_count() > 1 and multiple_gpu:
            self.multiple_gpu=True
            model = nn.SyncBatchNorm.convert_sync_batchnorm(model)
            model = nn.DataParallel(model)
        else:
            self.multiple_gpu = False
        self.device = device
        self.model = model.to(device)

        # model.to(device)
        self.data = data
        self.batch_size = batch_size
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        self.losses = []

    def train(self):
        # with tqdm(total=self.data.shape[0], desc="Training model", unit="%", unit_scale=True) as pbar:
        #     pbar.update(self.batch_size)
        train_loader = DataLoader(self.data, batch_size=self.batch_size, shuffle=True)
        for epoch in  tqdm(range(self.epochs),desc="Training model ...", ):
            epoch_loss = 0
            for x_batch, mask_batch in train_loader:
                x_batch, mask_batch = x_batch.to(self.device), mask_batch.to(self.device)
                self.optimizer.zero_grad()
                recon_x, mu, log_var = self.model(x_batch, mask_batch)
                if self.multiple_gpu:
                    loss = self.model.module.loss_function(recon_x, x_batch, mu, log_var, mask_batch)
                else:
                    loss = self.model.loss_function(recon_x, x_batch, mu, log_var, mask_batch)
                loss.backward()
                self.optimizer.step()
                epoch_loss += loss.item()
        avg_loss = epoch_loss / len(train_loader)
        self.losses.append(avg_loss)
        print(f"Epoch {epoch + 1}/{self.epochs}, Loss: {avg_loss:.4f}")
        if epoch>5 and avg_loss<0.0:
            return


    def save_model(self,path_model):
        if self.multiple_gpu:
            torch.save(self.model.module, f"{path_model}")
        else:
            torch.save(self.model, f"{path_model}")
    def get_latent_features(self, data):
        if self.multiple_gpu:
            return self.model.module.generate_latent_features(data.to(self.device))
        else:
            return self.model.generate_latent_features(data.to(self.device))

    def get_losses(self):
        return self.losses