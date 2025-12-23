import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# 定义目标网格的经度和纬度范围
lon_range = np.arange(-180, 180, 2)
lat_range = np.arange(-60, 92, 2)

# 定义重采样函数
def bilinear_interpolation(data):
    # 使用xarray的interp函数进行双线性插值
    resampled_data = data.interp(lon=lon_range, lat=lat_range, method='linear')
    return resampled_data


def adjust_longitude(dataset):
    lon_name = 'lon'
    dataset = dataset.assign_coords({lon_name: (((dataset[lon_name] + 180) % 360) - 180)})
    dataset = dataset.sortby(lon_name)
    
    # 检查是否已经存在 180 度经度值，如果没有则添加
    if 180 not in dataset[lon_name].values:
        # 找到 -180 度的数据
        minus_180_data = dataset.sel({lon_name: -180}, method="nearest")
        
        # 创建新的经度数组（包含 180 度）
        new_lon_values = np.append(dataset[lon_name].values, 180)
        
        # 再次排序，确保单调递增
        new_lon_values_sorted = np.sort(new_lon_values)
        
        # 扩展数据集（使用 reindex + fillna 确保数据连续性）
        dataset = dataset.reindex(
            {lon_name: new_lon_values_sorted},
            method="nearest",  # 使用最近邻填充
            tolerance=10  # 允许的最大距离（可选）
        )
        
        # 确保 180 度的数据与 -180 度相同
        dataset.loc[{lon_name: 180}] = minus_180_data
    
    return dataset


# # 用于UKESM模式（-179.1, 179.1）
# def adjust_longitude(dataset):
#     lon_name = 'lon'
#     dataset = dataset.assign_coords({lon_name: (((dataset[lon_name] + 180) % 360) - 180)})
#     dataset = dataset.sortby(lon_name)
    
#     # 检查是否已经包含 -180 和 180，如果没有则扩展
#     lon_min = dataset[lon_name].min().item()
#     lon_max = dataset[lon_name].max().item()
    
#     if lon_min > -180 or lon_max < 180:
#         # 获取边界数据（用于填充 -180 和 180）
#         left_data = dataset.sel({lon_name: lon_min}, method="nearest")
#         right_data = dataset.sel({lon_name: lon_max}, method="nearest")
        
#         # 创建新的经度数组（包含 -180 和 180）
#         new_lon_values = np.unique(
#             np.concatenate([
#                 [-180],  # 强制添加 -180
#                 dataset[lon_name].values,
#                 [180]    # 强制添加 180
#             ])
#         )
        
#         # 重新索引数据
#         dataset = dataset.reindex(
#             {lon_name: new_lon_values},
#             method="nearest",  
#             tolerance=10       # 允许的最大距离（可选）
#         )
        
#         # 确保 -180 和 180 的数据与边界数据一致
#         dataset.loc[{lon_name: -180}] = right_data
#         dataset.loc[{lon_name: 180}] = left_data
    
#     return dataset


f_spei = r"E:\SPEI_results\ACCESS\1all_SPEI3.nc"
spei_30 = xr.open_dataset(f_spei).spei.isel(time=slice(12, 372))
hou30 = xr.open_dataset(f_spei).spei.isel(time=slice(1260, 1620))
spei_tas1 = xr.open_dataset(f_spei).spei.isel(time=slice(1632, 1992))
spei_tas2 = xr.open_dataset(f_spei).spei.isel(time=slice(2004, 2364))
spei_tas3 = xr.open_dataset(f_spei).spei.isel(time=slice(2376, 2736))
spei_rn1 = xr.open_dataset(f_spei).spei.isel(time=slice(2748, 3108))
spei_rn2 = xr.open_dataset(f_spei).spei.isel(time=slice(3120, 3480))
spei_rn3 = xr.open_dataset(f_spei).spei.isel(time=slice(3492, 3852))
spei_rh1 = xr.open_dataset(f_spei).spei.isel(time=slice(3864, 4224))
spei_rh2 = xr.open_dataset(f_spei).spei.isel(time=slice(4236, 4596))
spei_rh3 = xr.open_dataset(f_spei).spei.isel(time=slice(4608, 4968))
spei_u21 = xr.open_dataset(f_spei).spei.isel(time=slice(4980, 5340))
spei_u22 = xr.open_dataset(f_spei).spei.isel(time=slice(5352, 5712))
spei_u23 = xr.open_dataset(f_spei).spei.isel(time=slice(5724, 6084))
spei_pr1 = xr.open_dataset(f_spei).spei.isel(time=slice(6096, 6456))
spei_pr2 = xr.open_dataset(f_spei).spei.isel(time=slice(6468, 6828))
spei_pr3 = xr.open_dataset(f_spei).spei.isel(time=slice(6840, 7200))
spei_co2 = xr.open_dataset(f_spei).spei.isel(time=slice(7212, 7572))

# 计算时间平均SPEI3
SPEI3 = spei_30.mean(dim='time')
Hou30 = hou30.mean(dim='time')
Stas1 = spei_tas1.mean(dim='time')
Stas2 = spei_tas2.mean(dim='time')
Stas3 = spei_tas3.mean(dim='time')
Srn1 = spei_rn1.mean(dim='time')
Srn2 = spei_rn2.mean(dim='time')
Srn3 = spei_rn3.mean(dim='time')
Srh1 = spei_rh1.mean(dim='time')
Srh2 = spei_rh2.mean(dim='time')
Srh3 = spei_rh3.mean(dim='time')
Su21 = spei_u21.mean(dim='time')
Su22 = spei_u22.mean(dim='time')
Su23 = spei_u23.mean(dim='time')
Spr1 = spei_pr1.mean(dim='time')
Spr2 = spei_pr2.mean(dim='time')
Spr3 = spei_pr3.mean(dim='time')
Sco2 = spei_co2.mean(dim='time')

# 将每个时间平均后的结果转换为一维数组
SPEI3_flat = SPEI3.values.flatten()
Hou30_flat = Hou30.values.flatten()
Stas1_flat = Stas1.values.flatten()
Stas2_flat = Stas2.values.flatten()
Stas3_flat = Stas3.values.flatten()
Srn1_flat = Srn1.values.flatten()
Srn2_flat = Srn2.values.flatten()
Srn3_flat = Srn3.values.flatten()
Srh1_flat = Srh1.values.flatten()
Srh2_flat = Srh2.values.flatten()
Srh3_flat = Srh3.values.flatten()
Su21_flat = Su21.values.flatten()
Su22_flat = Su22.values.flatten()
Su23_flat = Su23.values.flatten()
Spr1_flat = Spr1.values.flatten()
Spr2_flat = Spr2.values.flatten()
Spr3_flat = Spr3.values.flatten()
Sco2_flat = Sco2.values.flatten()

# 定义X1-11
X1 = Stas1_flat - SPEI3_flat 
X2 = Stas2_flat - SPEI3_flat 
X3 = Stas3_flat - SPEI3_flat 
X4 = Srn1_flat - SPEI3_flat
X5 = Srn2_flat - SPEI3_flat 
X6 = Srn3_flat - SPEI3_flat
X7 = Srh1_flat - SPEI3_flat
X8 = Srh2_flat - SPEI3_flat
X9 = Srh3_flat - SPEI3_flat
X10 = Su21_flat - SPEI3_flat
X11 = Su22_flat - SPEI3_flat
X12 = Su23_flat - SPEI3_flat
X13 = Spr1_flat - SPEI3_flat
X14 = Spr2_flat - SPEI3_flat
X15 = Spr3_flat - SPEI3_flat
X16 = Sco2_flat - SPEI3_flat

# 构造方程组的系数矩阵
A = np.array([[0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]])

# 构造方程右侧的常数向量
X = np.array([X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, X16])

# 求解方程组
result = np.linalg.solve(A, X)
C_T1 = result[0]  
C_T2 = result[1]  
C_T3 = result[2] 
C_RN1 = result[3] 
C_RN2 = result[4]  
C_RN3 = result[5] 
C_RH1 = result[6]  
C_RH2 = result[7] 
C_RH3 = result[8] 
C_U1 = result[9]  
C_U2 = result[10] 
C_U3 = result[11] 
C_PRE1 = result[12] 
C_PRE2 = result[13] 
C_PRE3 = result[14]
C_CO2 = result[15]
C_all = Hou30_flat - SPEI3_flat


# 创建结果数据集并保存为NetCDF文件
def save_as_netcdf(data, var_name, filename):
    lon = SPEI3.lon
    lat = SPEI3.lat
    
    ds = xr.Dataset(
        {var_name: (('lat', 'lon'), data.reshape(len(lat), len(lon)))},
        coords={'lat': lat, 'lon': lon}
    )
    
    data_ad = adjust_longitude(ds)
    rdata = bilinear_interpolation(data_ad)
    rdata.to_netcdf(filename)

# 保存各个分量结果
save_as_netcdf(C_T1, 'T1', 'E:\\SPEI_results\\ACCESS\\1C_T1.nc')
save_as_netcdf(C_T2, 'T2', 'E:\\SPEI_results\\ACCESS\\1C_T2.nc')
save_as_netcdf(C_T3, 'T3', 'E:\\SPEI_results\\ACCESS\\1C_T3.nc')
save_as_netcdf(C_RN1, 'RN1', 'E:\\SPEI_results\\ACCESS\\1C_RN1.nc')
save_as_netcdf(C_RN2, 'RN2', 'E:\\SPEI_results\\ACCESS\\1C_RN2.nc')
save_as_netcdf(C_RN3, 'RN3', 'E:\\SPEI_results\\ACCESS\\1C_RN3.nc')
save_as_netcdf(C_RH1, 'RH1', 'E:\\SPEI_results\\ACCESS\\1C_RH1.nc')
save_as_netcdf(C_RH2, 'RH2', 'E:\\SPEI_results\\ACCESS\\1C_RH2.nc')
save_as_netcdf(C_RH3, 'RH3', 'E:\\SPEI_results\\ACCESS\\1C_RH3.nc')
save_as_netcdf(C_U1, 'U1', 'E:\\SPEI_results\\ACCESS\\1C_U1.nc')
save_as_netcdf(C_U2, 'U2', 'E:\\SPEI_results\\ACCESS\\1C_U2.nc')
save_as_netcdf(C_U3, 'U3', 'E:\\SPEI_results\\ACCESS\\1C_U3.nc')
save_as_netcdf(C_PRE1, 'PRE1', 'E:\\SPEI_results\\ACCESS\\1C_PRE1.nc')
save_as_netcdf(C_PRE2, 'PRE2', 'E:\\SPEI_results\\ACCESS\\1C_PRE2.nc')
save_as_netcdf(C_PRE3, 'PRE3', 'E:\\SPEI_results\\ACCESS\\1C_PRE3.nc')
save_as_netcdf(C_CO2, 'CO2', 'E:\\SPEI_results\\ACCESS\\1C_CO2.nc')
save_as_netcdf(C_all, 'CTL', 'E:\\SPEI_results\\ACCESS\\1C_all.nc')

