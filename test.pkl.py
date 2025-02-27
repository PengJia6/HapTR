import pickle
import argparse
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from torchvision import datasets, transforms
from torch import nn
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence


# 1. 定义VAE模型
class VAE(nn.Module):
    def __init__(self, input_dim, latent_dim, hidden_dim):
        super(VAE, self).__init__()
        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.hidden_dim = hidden_dim

        # 编码器：LSTM
        self.lstm = nn.LSTM(input_dim, hidden_dim, batch_first=True, bidirectional=True)

        # 通过LSTM的输出获得潜在空间的均值和对数方差
        self.fc_mu = nn.Linear(hidden_dim * 2, latent_dim)  # 由于是双向LSTM，所以乘2
        self.fc_logvar = nn.Linear(hidden_dim * 2, latent_dim)

        # 解码器：使用LSTM将潜在空间映射回原始序列的空间
        self.fc_decode = nn.Linear(latent_dim, hidden_dim * 2)
        self.lstm_decode = nn.LSTM(hidden_dim, input_dim, batch_first=True)

    def encode(self, x, lengths):
        # 对输入序列进行填充和打包
        # print(x,lengths)
        packed_input = pack_padded_sequence(x, lengths, batch_first=True, enforce_sorted=False)
        _, (hn, _) = self.lstm(packed_input)

        # 获取均值和对数方差
        mu = self.fc_mu(hn[-1])  # 使用最后一个隐状态
        logvar = self.fc_logvar(hn[-1])
        return mu, logvar

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, seq_length):
        # 解码器的输入
        z = self.fc_decode(z).unsqueeze(0).repeat(seq_length, 1, 1)  # 扩展维度
        output, _ = self.lstm_decode(z)
        return output

    def forward(self, x, lengths):
        mu, logvar = self.encode(x, lengths)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decode(z, x.size(1))  # 重建序列
        return recon_x, mu, logvar


# 2. 损失函数：重构误差 + KL散度
def vae_loss(recon_x, x, mu, logvar):
    # 使用MSE作为重建误差
    reconstruction_loss = nn.MSELoss(reduction='sum')(recon_x, x)

    # 计算KL散度
    kl_divergence = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    return reconstruction_loss + kl_divergence


# 3. 定义数据集类
class SequenceDataset(Dataset):
    def __init__(self, sequences, max_len):
        self.sequences = sequences
        self.max_len = max_len

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        # 填充至最大长度
        padded_sequence = torch.tensor(sequence + [[0, 0, 0]] * (self.max_len - len(sequence)), dtype=torch.float32)
        return padded_sequence


# 4. 训练VAE模型
def train_vae(dataloader, vae, optimizer, num_epochs=10):
    vae.train()
    for epoch in range(num_epochs):
        running_loss = 0.0
        for data in dataloader:
            # 获取数据和长度信息
            sequences = data
            # lengths = [i[0] for i in torch.sum(sequences != 0, dim=1)]  # 计算每个序列的实际长度
            lengths = torch.sum(sequences != 0, dim=1)  # 计算每个序列的实际长度
            print(lengths)
            print(sequences)
            optimizer.zero_grad()
            recon_batch, mu, logvar = vae(sequences, lengths)
            loss = vae_loss(recon_batch, sequences, mu, logvar)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()

        print(f"Epoch {epoch + 1}/{num_epochs}, Loss: {running_loss / len(dataloader)}")


# 8. 特征提取
def extract_features(vae, dataloader):
    vae.eval()
    features = []
    with torch.no_grad():
        for data in dataloader:
            lengths = torch.sum(data != 0, dim=1)
            mu, _ = vae.encode(data, lengths)
            features.append(mu.cpu().numpy())
    return np.concatenate(features, axis=0)


def train():
    # 5. 示例数据生成
    sequences = [torch.rand(np.random.randint(5, 20), 3).tolist() for _ in range(100)]  # 100个序列，维度为3的多维序列
    # print([for i in seq])
    max_len = max([len(seq) for seq in sequences])  # 找到最大长度
    # print(max_len)
    # 6. 数据加载
    dataset = SequenceDataset(sequences, max_len)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

    # 7. 模型训练
    vae = VAE(input_dim=3, latent_dim=5, hidden_dim=16)
    optimizer = optim.Adam(vae.parameters(), lr=0.001)

    train_vae(dataloader, vae, optimizer)
    features = extract_features(vae, dataloader)
    print("Extracted features shape:", features.shape)




def calculate_similarity_matrix(seq1, seq2, window_size=3, match=2, mismatch=-1, gap=-2):
    """
    滑动窗口 + 动态规划算法，计算局部相似性矩阵
    :param seq1: DNA序列1
    :param seq2: DNA序列2
    :param window_size: 滑动窗口大小
    :param match: 匹配得分
    :param mismatch: 不匹配惩罚
    :param gap: 空位惩罚
    :return: 局部相似性矩阵
    """
    len1, len2 = len(seq1), len(seq2)
    similarity_matrix = np.zeros((len1 - window_size + 1, len2 - window_size + 1))

    # 定义一个局部比对的动态规划函数
    def local_alignment(subseq1, subseq2):
        m, n = len(subseq1), len(subseq2)
        dp = np.zeros((m + 1, n + 1))
        max_score = 0

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match_score = match if subseq1[i - 1] == subseq2[j - 1] else mismatch
                dp[i][j] = max(0,  # 局部比对，允许非负
                               dp[i - 1][j - 1] + match_score,  # 匹配/不匹配
                               dp[i - 1][j] + gap,  # 空位1
                               dp[i][j - 1] + gap)  # 空位2
                max_score = max(max_score, dp[i][j])
        return max_score

    # 滑动窗口遍历两个序列
    for i in range(len1 - window_size + 1):
        subseq1 = seq1[i:i + window_size]
        for j in range(len2 - window_size + 1):
            subseq2 = seq2[j:j + window_size]
            similarity_matrix[i][j] = local_alignment(subseq1, subseq2)

    return similarity_matrix


def read_pickle(input_file):
    with open(input_file, 'rb') as file:
        num = 1
        while True:
            try:
                num += 1
                region = pickle.load(file)
                print("读取的数据:", num)
                region_variant_info = region.variants_info
                region_repeat_info = region.repeats
                for idx, repeat in region_repeat_info.items():
                    # print(idx,repeat)
                    print(repeat.content, "===========================")
                    ref_similarity_matrix = calculate_similarity_matrix(repeat.content, repeat.content)
                    ref_matrix_mean = np.nanmean(ref_similarity_matrix, axis=0)
                    for read_id, read_features in repeat.repeat_feature.items():
                        if (read_features.seq_list is not None) and (read_features.mut_list is not None) and (read_features.qual_list is not None):
                            read_similarity_matrix = np.nanmean(calculate_similarity_matrix("".join(read_features.seq_list), repeat.content), axis=0) - ref_matrix_mean
                            print(read_similarity_matrix.shape, len(read_features.mut_list), len(read_features.seq_list), len(read_features.qual_list))

                            # print(read_id, len(read_features.mut_list), len(read_features.seq_list), len(read_features.qual_list))
                        # print(read_id,)
                        # print(read_id, read_features.qual_list)
            except EOFError:
                break


def main():
    parser = argparse.ArgumentParser(description="Extract the repeat sequence from assemblies.")
    parser.add_argument("-i", '--input', type=str, required=True, help="The input pickle file.")
    parser.add_argument("-o", '--output', required=True, type=str, help="The output of the benchmark. [Required]")
    # parser.add_argument('-r', "--repeat", required=True, type=str,
    #                     help="The path of tandem repeat file program. [Required]")
    # parser.add_argument("-b1", '--bam_hap1', type=str, required=True,
    #                     help="The aligned bam file of the hap1 fasta. [Required]")
    # parser.add_argument("-b2", '--bam_hap2', type=str, required=True,
    #                     help="The aligned bam file of the hap1 fasta. [Required]")
    # see the help documentation at https://tandem.bu.edu/trf/help
    args = parser.parse_args()
    read_pickle(args.input)


if __name__ == "__main__":
    main()
    # train()
