
import pickle

class xxx:
    def __init__(self,a):
        self.a=1
        self.b=a
# 要保存的数据
data = {"name": "Yuting Shang", "age": 30, "skills": ["Python", "Genomics", "Bioinformatics"]}
data2 = {"name2": "Yuting Shang", "age": 90, "skills": ["Python", "Genomics", "Bioinformatics"]}

# 写入 pkl 文件
with open("data.pkl", "wb") as file:  # 使用 'wb' 模式打开文件
    pickle.dump(xxx(data), file)
    pickle.dump(xxx(data2), file)

print("数据已保存到 data.pkl")
# 读取 pkl 文件
with open("data.pkl", "rb") as file:  # 使用 'rb' 模式打开文件
    while True:
        try:
            loaded_data = pickle.load(file)
            print(loaded_data.b,"lll")
        except EOFError:
            break
        # loaded_data2 = pickle.load(file)

# print("从 data.pkl 加载的数据：", loaded_data)
# print("从 data.pkl 加载的数据：", loaded_data2)
