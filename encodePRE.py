
def convert_to_integer_3d(point_cloud, target_max=2**15):
    # 将输入转换为 NumPy 数组
    point_cloud = np.array(point_cloud)
    
    # 获取数据的最小值和最大值
    min_val = point_cloud.min()
    max_val = point_cloud.max()
    
    # 统一归一化到 [0, 1]
    normalized = (point_cloud - min_val) / (max_val - min_val)
    
    # 缩放到目标整数范围
    scaled = normalized * (target_max - 1)  # 缩放到 [0, target_max)
    
    # 四舍五入并转换为整数
    int_data = np.round(scaled).astype(int)
    
    # 确保在目标范围内
    int_data = np.clip(int_data, 0, target_max - 1)  # 限定在 [0, target_max - 1] 的范围内
    
    return int_data

def read_point_cloud_from_file(filename):
    point_cloud = []
    with open(filename, 'r') as file:
        for line in file:
            # 将每行的字符串转换为浮点数并添加到列表中
            point = list(map(float, line.split()))
            point_cloud.append(point)
    return point_cloud

def write_point_cloud_to_file(filename, point_cloud):
    with open(filename, 'w') as file:
        for point in point_cloud:
            # 将每个点格式化为字符串并写入文件
            file.write(f"{point[0]} {point[1]} {point[2]}\n")








def encode_value(coordinate, x):
    ans = 0
    position = 0
    for idx in x:
        idx = int(idx)
       # Debug: Print the current index and coordinate value
       # print(f"Index: {idx}, Coordinate Value: {coordinate[idx]}, Position: {position}, Ans before: {ans}")
        ans |= (coordinate[idx] & 1) << position
        coordinate[idx] >>= 1
       # Debug: Print the updated answer
       # print(f"Ans after: {ans}")
        position += 1
    return ans

def read_coordinates_from_file(filename):
    coordinates = []
    with open(filename, 'r') as file:
        for line in file:
            values = list(map(int, line.split()))
            coordinates.extend(values)
    return coordinates

def write_encoded_values_to_file(filename, encoded_values):
    with open(filename, 'w') as file:
        for value in encoded_values:
            file.write(f"{value}\n")










import numpy as np

class Point:
    def __init__(self, x, y, z, v=0.0):
        self.x = x
        self.y = y
        self.z = z
        self.v = v

def read_points_from_file_a(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            values = list(map(float, line.split()))
            if len(values) >= 3:
                point = Point(values[0], values[1], values[2])
                points.append(point)
    return points

def read_v_values_from_file_b(filename, points):
    with open(filename, 'r') as file:
        for index, line in enumerate(file):
            if index < len(points):
                points[index].v = float(line.strip())

def write_sorted_points_to_file(filename, points):
    with open(filename, 'w') as file:
        for point in points:
            file.write(f"{point.x:.15f} {point.y:.15f} {point.z:.15f}\n")

def encoder(path):
    # 示例使用
    input_filename = 'PRE/pre_8192/' +path  # 假设你的点云数据保存在这个文件中
    output_filename = 'int_points_random.txt'  # 输出文件名
    print(input_filename)
# 读取点云数据
    point_cloud = read_point_cloud_from_file(input_filename)

# 转换为整数点云数据
    int_point_cloud = convert_to_integer_3d(point_cloud)

# 将结果写入文件
    write_point_cloud_to_file(output_filename, int_point_cloud)

    print(f"Converted integer point cloud saved to {output_filename}")
    encoded_values = []
    with open("int_points_random.txt", 'r') as input_file:
        for line in input_file:
            coordinates = list(map(int, line.split()))
            
            # Fixed x vector
            x = [0, 1, 2] * 15

            # Debug: Print the coordinates being processed
            # print(f"Coordinates: {coordinates}")

            # Encode current line of coordinates
            ans = encode_value(coordinates, x)
            encoded_values.append(ans)

    # Write encoded values to file
    write_encoded_values_to_file("encoded_values.txt", encoded_values)
    
    print("Encoded values saved to encoded_values.txt")

        # Read data from files
    points = read_points_from_file_a(input_filename)
    read_v_values_from_file_b("encoded_values.txt", points)

    # Sort points by v value
    points.sort(key=lambda point: point.v)

    # Write sorted results to file
    write_sorted_points_to_file('NEWPRE/pre_8192/'+path, points)
    
    print("Sorted points saved to sorted_points.txt")


    # 定义读取文件的函数
# 定义读取文件的函数
# 定义读取文件的函数
def read_file_lines_with_index(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                for i in range(0, 8):  # i 从 0 到 8
                    path = line.strip() +  "_" + str(i) + ".txt"  # 添加索引和".txt"
                    encoder(path)
                    # print(path)
    except FileNotFoundError:
        print(f"文件 {file_path} 未找到。")
    except Exception as e:
        print(f"发生错误: {e}")


# 替换为你的文件路径
read_file_lines_with_index("compare.list")