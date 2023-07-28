import torch

def string_to_int_list(s, max_len=50):
    int_list = [ord(c) for c in s[:max_len]]
    if len(int_list) < max_len:
        int_list += [0] * (max_len - len(int_list))
    return int_list

# 示例
strings = ["hello1212", "world"]
max_len = max(len(s) for s in strings)

tensor = torch.tensor([string_to_int_list(s) for s in strings], dtype=torch.int8)

print(tensor)