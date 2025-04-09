# ======================================================================================================================
# Project: HapTR
# Script : test.smk TODO check 
# Author : Peng Jia
# Date   :  2025/3/4
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
print("ll")
with open("/data/home/pengjia/HapTR/train_test/LCL5.HiFi.v0.1.4_test.model.importance.yaml") as f:
    data=yaml.safe_load(f)
print(data)
x=[]
y=[]
eval_num=data["eval_num"]
data_importance=data["importance"]
for i in range(200):
    x.append(i)

    if i in data_importance:
        y.append(data_importance[i]/eval_num)
    else:
        y.append(0)
plt.bar(x,y)
plt.show()
sns.kdeplot(y)
plt.show()
# print(len([i for i in y if i >0.2]))
print()
print(data_importance)
