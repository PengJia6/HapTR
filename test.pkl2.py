import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from torch.utils.data import DataLoader, TensorDataset
from torch.nn.utils.rnn import pad_sequence
import torch.nn.functional as F

# 设置设备



class VAEModel(nn.Module):
    def __init__(self, timesteps, input_dim, latent_dim=2, lstm_units=128):
        super(VAEModel, self).__init__()
        self.timesteps = timesteps
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.lstm_units = lstm_units

        # 编码器部分
        self.encoder_lstm = nn.LSTM(input_dim, lstm_units, batch_first=True)
        self.fc_mu = nn.Linear(lstm_units, latent_dim)
        self.fc_log_var = nn.Linear(lstm_units, latent_dim)

        # 解码器部分
        self.decoder_lstm = nn.LSTM(latent_dim, lstm_units, batch_first=True)
        self.decoder_fc = nn.Linear(lstm_units, input_dim)

    def reparameterize(self, mu, log_var):
        # 重参数化技巧
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x, mask):
        # 编码器
        lstm_out, _ = self.encoder_lstm(x)
        mu = self.fc_mu(lstm_out[:, -1, :])  # 取最后一个时间步的输出
        log_var = self.fc_log_var(lstm_out[:, -1, :])

        # 重参数化
        z = self.reparameterize(mu, log_var)

        # 解码器
        z = z.unsqueeze(1).repeat(1, self.timesteps, 1)  # 扩展到timesteps
        lstm_out, _ = self.decoder_lstm(z)
        recon_x = self.decoder_fc(lstm_out)

        # 只计算mask部分的重构损失
        recon_x = recon_x * mask.unsqueeze(-1)  # 将mask应用于输出

        return recon_x, mu, log_var

    def loss_function(self, recon_x, x, mu, log_var, mask):
        # 计算重构损失（MSE）+ KL散度损失
        reconstruction_loss = F.mse_loss(recon_x, x, reduction='none')  # 不进行汇总，保留每个时间步的损失
        reconstruction_loss = torch.sum(reconstruction_loss * mask.unsqueeze(-1))  # 只计算非填充部分的损失
        kl_divergence = -0.5 * torch.mean(1 + log_var - mu.pow(2) - log_var.exp())
        return reconstruction_loss + kl_divergence

    def generate_latent_features(self, x):
        # 生成潜在空间的特征
        lstm_out, _ = self.encoder_lstm(x)
        mu = self.fc_mu(lstm_out[:, -1, :])  # 取最后一个时间步的输出
        return mu


class VAETrainer:
    def __init__(self, model, data, batch_size=32, epochs=50, learning_rate=1e-3):
        self.model = model.to(device)  # 将模型移至GPU
        self.data = data
        self.batch_size = batch_size
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        self.losses = []

    def train(self):
        # 创建数据加载器
        train_loader = DataLoader(self.data, batch_size=self.batch_size, shuffle=True)

        for epoch in range(self.epochs):
            epoch_loss = 0
            for x_batch, mask_batch in train_loader:
                x_batch, mask_batch = x_batch.to(device), mask_batch.to(device)  # 将数据移至GPU
                self.optimizer.zero_grad()

                # 计算模型输出
                recon_x, mu, log_var = self.model(x_batch, mask_batch)

                # 计算损失
                loss = self.model.loss_function(recon_x, x_batch, mu, log_var, mask_batch)
                loss.backward()

                # 更新参数
                self.optimizer.step()

                epoch_loss += loss.item()

            average_loss = epoch_loss / len(train_loader)
            self.losses.append(average_loss)
            print(f"Epoch {epoch + 1}/{self.epochs}, Loss: {average_loss:.4f}")

    def get_latent_features(self, data):
        # 提取潜在特征
        return self.model.generate_latent_features(data.to(device))  # 将数据移至GPU

    def get_losses(self):
        return self.losses

    @staticmethod
    def generate_simulated_data(num_sequences, timesteps, input_dim):
        # 生成模拟数据
        sequences = np.random.rand(num_sequences, timesteps, input_dim)
        lengths = np.random.randint(1, timesteps + 1, num_sequences)  # 随机生成不等长序列
        mask = np.zeros((num_sequences, timesteps), dtype=np.float32)
        for i, length in enumerate(lengths):
            mask[i, :length] = 1  # 为有效部分设置mask
        return torch.tensor(sequences, dtype=torch.float32), torch.tensor(mask, dtype=torch.float32)


# 设置随机种子以保证结果可复现
np.random.seed(42)
torch.manual_seed(42)

# 模拟数据参数
num_sequences = 1000  # 生成100个序列
timesteps = 5  # 每个序列5个时间步
input_dim = 3  # 每个时间步3个特征

# 生成模拟数据
sequences, masks = VAETrainer.generate_simulated_data(num_sequences, timesteps, input_dim)
print("sequences", sequences.shape)
print("mask", masks.shape)
print(sequences)
# 创建VAE模型实例
vae_model = VAEModel(timesteps=timesteps, input_dim=input_dim)

# 创建训练器实例
trainer = VAETrainer(model=vae_model, data=TensorDataset(sequences, masks), epochs=50, batch_size=32)

# 训练VAE模型
trainer.train()

# 提取潜在特征
latent_features = trainer.get_latent_features(sequences)

# 打印潜在特征的形状
print("Latent Features Shape:", latent_features.shape)

# 获取训练过程中每个epoch的损失
losses = trainer.get_losses()
print("Losses over epochs:", losses)
