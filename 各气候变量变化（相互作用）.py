import numpy as np
import xarray as xr
from math import pi, sin, cos, tan, log

# 定义目标网格的经度和纬度范围
lon_range = np.arange(-180, 180, 2)
lat_range = np.arange(-60, 92, 2)

# 定义重采样函数
def bilinear_interpolation(data):
    # 使用xarray的interp函数进行双线性插值
    resampled_data = data.interp(lon=lon_range, lat=lat_range, method='linear')
    return resampled_data


# def adjust_longitude(dataset):
#     lon_name = 'lon'
#     dataset = dataset.assign_coords({lon_name: (((dataset[lon_name] + 180) % 360) - 180)})
#     dataset = dataset.sortby(lon_name)
    
#     # 检查是否已经存在 180 度经度值，如果没有则添加
#     if 180 not in dataset[lon_name].values:
#         # 找到 -180 度的数据
#         minus_180_data = dataset.sel({lon_name: -180}, method="nearest")
        
#         # 创建新的经度数组（包含 180 度）
#         new_lon_values = np.append(dataset[lon_name].values, 180)
        
#         # 再次排序，确保单调递增
#         new_lon_values_sorted = np.sort(new_lon_values)
        
#         # 扩展数据集（使用 reindex + fillna 确保数据连续性）
#         dataset = dataset.reindex(
#             {lon_name: new_lon_values_sorted},
#             method="nearest",  # 使用最近邻填充
#             tolerance=10  # 允许的最大距离（可选）
#         )
        
#         # 确保 180 度的数据与 -180 度相同
#         dataset.loc[{lon_name: 180}] = minus_180_data
    
#     return dataset


# 用于UKESM模式（-179.1, 179.1）
def adjust_longitude(dataset):
    lon_name = 'lon'
    dataset = dataset.assign_coords({lon_name: (((dataset[lon_name] + 180) % 360) - 180)})
    dataset = dataset.sortby(lon_name)
    
    # 检查是否已经包含 -180 和 180，如果没有则扩展
    lon_min = dataset[lon_name].min().item()
    lon_max = dataset[lon_name].max().item()
    
    if lon_min > -180 or lon_max < 180:
        # 获取边界数据（用于填充 -180 和 180）
        left_data = dataset.sel({lon_name: lon_min}, method="nearest")
        right_data = dataset.sel({lon_name: lon_max}, method="nearest")
        
        # 创建新的经度数组（包含 -180 和 180）
        new_lon_values = np.unique(
            np.concatenate([
                [-180],  # 强制添加 -180
                dataset[lon_name].values,
                [180]    # 强制添加 180
            ])
        )
        
        # 重新索引数据
        dataset = dataset.reindex(
            {lon_name: new_lon_values},
            method="nearest",  
            tolerance=10       # 允许的最大距离（可选）
        )
        
        # 确保 -180 和 180 的数据与边界数据一致
        dataset.loc[{lon_name: -180}] = right_data
        dataset.loc[{lon_name: 180}] = left_data
    
    return dataset


# ##### 降水变化
# f_prCTL = r"E:\CMIP6\UKESM\pr_1pctCO2_1850-1999.nc"
# prCTLq30 = xr.open_dataset(f_prCTL).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
# prCTLh30 = xr.open_dataset(f_prCTL).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24
# f_prrad = r"E:\CMIP6\UKESM\pr_1pctCO2-rad_1850-1999.nc"
# prradq30 = xr.open_dataset(f_prrad).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
# prradh30 = xr.open_dataset(f_prrad).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24
# f_prbgc = r"E:\CMIP6\UKESM\pr_1pctCO2-bgc_1850-1999.nc"
# prbgcq30 = xr.open_dataset(f_prbgc).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
# prbgch30 = xr.open_dataset(f_prbgc).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24

# month_q30 = prCTLq30.time.dt.days_in_month
# month_h30 = prCTLh30.time.dt.days_in_month

# prCTLq30_total = prCTLq30 * month_q30
# prCTLh30_total = prCTLh30 * month_h30
# prradq30_total = prradq30 * month_q30
# prradh30_total = prradh30 * month_h30
# prbgcq30_total = prbgcq30 * month_q30
# prbgch30_total = prbgch30 * month_h30

# PCTL = prCTLh30_total.groupby('time.year').sum(dim='time').mean(dim='year') - prCTLq30_total.groupby('time.year').sum(dim='time').mean(dim='year')
# Prad = prradh30_total.groupby('time.year').sum(dim='time').mean(dim='year') - prradq30_total.groupby('time.year').sum(dim='time').mean(dim='year')
# Pbgc = prbgch30_total.groupby('time.year').sum(dim='time').mean(dim='year') - prbgcq30_total.groupby('time.year').sum(dim='time').mean(dim='year')

# xianghu = PCTL - Prad - Pbgc



# ##### 温度变化
# f_TmaxCTL = r"E:\CMIP6\ACCESS\tasmax_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
# TmaxCTLq30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('0102-01', '0131-12')) - 273.16
# TmaxCTLh30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('0206-01', '0235-12')) - 273.16
# f_TminCTL = r"E:\CMIP6\ACCESS\tasmin_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
# TminCTLq30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('0102-01', '0131-12')) - 273.16
# TminCTLh30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('0206-01', '0235-12')) - 273.16

# TCTLq30 = (TmaxCTLq30 + TminCTLq30)/2
# TCTLh30 = (TmaxCTLh30 + TminCTLh30)/2

# f_Tmaxrad = r"E:\CMIP6\ACCESS\tasmax_Amon_ACCESS-ESM1-5_1pctCO2-rad_r1i1p1f1_gn_010101-025012.nc"
# Tmaxradq30 = xr.open_dataset(f_Tmaxrad).tasmax.sel(time=slice('0102-01', '0131-12')) - 273.16
# Tmaxradh30 = xr.open_dataset(f_Tmaxrad).tasmax.sel(time=slice('0206-01', '0235-12')) - 273.16
# f_Tminrad = r"E:\CMIP6\ACCESS\tasmin_Amon_ACCESS-ESM1-5_1pctCO2-rad_r1i1p1f1_gn_010101-025012.nc"
# Tminradq30 = xr.open_dataset(f_Tminrad).tasmin.sel(time=slice('0102-01', '0131-12')) - 273.16
# Tminradh30 = xr.open_dataset(f_Tminrad).tasmin.sel(time=slice('0206-01', '0235-12')) - 273.16

# Tradq30 = (Tmaxradq30 + Tminradq30)/2
# Tradh30 = (Tmaxradh30 + Tminradh30)/2

# f_Tmaxbgc = r"E:\CMIP6\ACCESS\tasmax_Amon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
# Tmaxbgcq30 = xr.open_dataset(f_Tmaxbgc).tasmax.sel(time=slice('0102-01', '0131-12')) - 273.16
# Tmaxbgch30 = xr.open_dataset(f_Tmaxbgc).tasmax.sel(time=slice('0206-01', '0235-12')) - 273.16
# f_Tminbgc = r"E:\CMIP6\ACCESS\tasmin_Amon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
# Tminbgcq30 = xr.open_dataset(f_Tminbgc).tasmin.sel(time=slice('0102-01', '0131-12')) - 273.16
# Tminbgch30 = xr.open_dataset(f_Tminbgc).tasmin.sel(time=slice('0206-01', '0235-12')) - 273.16

# Tbgcq30 = (Tmaxbgcq30 + Tminbgcq30)/2
# Tbgch30 = (Tmaxbgch30 + Tminbgch30)/2

# TCTL = TCTLh30.mean(dim='time') - TCTLq30.mean(dim='time')
# Trad = Tradh30.mean(dim='time') - Tradq30.mean(dim='time')
# Tbgc = Tbgch30.mean(dim='time') - Tbgcq30.mean(dim='time')

# xianghu = TCTL - Trad - Tbgc



##### 净辐射变化
f_rssCTL = r"E:\CMIP6\UKESM\rss_1pctCO2_1850-1999.nc"
RnsCTLq30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1851-01', '1880-12'))
RnsCTLh30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1955-01', '1984-12'))
f_rlsCTL = r"E:\CMIP6\UKESM\rls_1pctCO2_1850-1999.nc"
RnlCTLq30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1851-01', '1880-12'))
RnlCTLh30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1955-01', '1984-12'))

month_q30 = RnsCTLq30.time.dt.days_in_month
month_h30 = RnsCTLh30.time.dt.days_in_month

RnCTLq30 = (RnsCTLq30 + RnlCTLq30) * 0.0864 * month_q30
RnCTLh30 = (RnsCTLh30 + RnlCTLh30) * 0.0864 * month_h30
RNCTLq30 = xr.where(RnCTLq30 < 0, 0, RnCTLq30)
RNCTLh30 = xr.where(RnCTLh30 < 0, 0, RnCTLh30)

f_rssrad = r"E:\CMIP6\UKESM\rss_1pctCO2-rad_1850-1999.nc"
Rnsradq30 = xr.open_dataset(f_rssrad).rss.sel(time=slice('1851-01', '1880-12'))
Rnsradh30 = xr.open_dataset(f_rssrad).rss.sel(time=slice('1955-01', '1984-12'))
f_rlsrad = r"E:\CMIP6\UKESM\rls_1pctCO2-rad_1850-1999.nc"
Rnlradq30 = xr.open_dataset(f_rlsrad).rls.sel(time=slice('1851-01', '1880-12'))
Rnlradh30 = xr.open_dataset(f_rlsrad).rls.sel(time=slice('1955-01', '1984-12'))

Rnradq30 = (Rnsradq30 + Rnlradq30) * 0.0864 * month_q30
Rnradh30 = (Rnsradh30 + Rnlradh30) * 0.0864 * month_h30
RNradq30 = xr.where(Rnradq30 < 0, 0, Rnradq30)
RNradh30 = xr.where(Rnradh30 < 0, 0, Rnradh30)

f_rssbgc = r"E:\CMIP6\UKESM\rss_1pctCO2-bgc_1850-1999.nc"
Rnsbgcq30 = xr.open_dataset(f_rssbgc).rss.sel(time=slice('1851-01', '1880-12'))
Rnsbgch30 = xr.open_dataset(f_rssbgc).rss.sel(time=slice('1955-01', '1984-12'))
f_rlsbgc = r"E:\CMIP6\UKESM\rls_1pctCO2-bgc_1850-1999.nc"
Rnlbgcq30 = xr.open_dataset(f_rlsbgc).rls.sel(time=slice('1851-01', '1880-12'))
Rnlbgch30 = xr.open_dataset(f_rlsbgc).rls.sel(time=slice('1955-01', '1984-12'))

Rnbgcq30 = (Rnsbgcq30 + Rnlbgcq30) * 0.0864 * month_q30
Rnbgch30 = (Rnsbgch30 + Rnlbgch30) * 0.0864 * month_h30
RNbgcq30 = xr.where(Rnbgcq30 < 0, 0, Rnbgcq30)
RNbgch30 = xr.where(Rnbgch30 < 0, 0, Rnbgch30)

RNCTL = RNCTLh30.groupby('time.year').sum(dim='time').mean(dim='year') - RNCTLq30.groupby('time.year').sum(dim='time').mean(dim='year')
RNrad = RNradh30.groupby('time.year').sum(dim='time').mean(dim='year') - RNradq30.groupby('time.year').sum(dim='time').mean(dim='year')
RNbgc = RNbgch30.groupby('time.year').sum(dim='time').mean(dim='year') - RNbgcq30.groupby('time.year').sum(dim='time').mean(dim='year')

xianghu = RNCTL - RNrad - RNbgc



# ##### 相对湿度变化
# f_RHCTL = r"E:\CMIP6\UKESM\hurs_1pctCO2_1850-1999.nc"
# RHCTLq30 = xr.open_dataset(f_RHCTL).hurs.sel(time=slice('1851-01', '1880-12'))
# RHCTLh30 = xr.open_dataset(f_RHCTL).hurs.sel(time=slice('1955-01', '1984-12'))
# f_RHrad = r"E:\CMIP6\UKESM\hurs_1pctCO2-rad_1850-1999.nc"
# RHradq30 = xr.open_dataset(f_RHrad).hurs.sel(time=slice('1851-01', '1880-12'))
# RHradh30 = xr.open_dataset(f_RHrad).hurs.sel(time=slice('1955-01', '1984-12'))
# f_RHbgc = r"E:\CMIP6\UKESM\hurs_1pctCO2-bgc_1850-1999.nc"
# RHbgcq30 = xr.open_dataset(f_RHbgc).hurs.sel(time=slice('1851-01', '1880-12'))
# RHbgch30 = xr.open_dataset(f_RHbgc).hurs.sel(time=slice('1955-01', '1984-12'))

# RHCTL = RHCTLh30.mean(dim='time') - RHCTLq30.mean(dim='time')
# RHrad = RHradh30.mean(dim='time') - RHradq30.mean(dim='time')
# RHbgc = RHbgch30.mean(dim='time') - RHbgcq30.mean(dim='time')

# xianghu = RHCTL - RHrad - RHbgc



# ##### 风速变化
# f_UCTL = r"E:\CMIP6\ACCESS\sfcWind_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc"
# UCTLq30 = xr.open_dataset(f_UCTL).sfcWind.sel(time=slice('0102-01', '0131-12')) * (4.87 / log(67.8 * 10 - 5.42))
# UCTLh30 = xr.open_dataset(f_UCTL).sfcWind.sel(time=slice('0206-01', '0235-12')) * (4.87 / log(67.8 * 10 - 5.42))
# f_Urad = r"E:\CMIP6\ACCESS\sfcWind_Amon_ACCESS-ESM1-5_1pctCO2-rad_r1i1p1f1_gn_010101-025012.nc"
# Uradq30 = xr.open_dataset(f_Urad).sfcWind.sel(time=slice('0102-01', '0131-12')) * (4.87 / log(67.8 * 10 - 5.42))
# Uradh30 = xr.open_dataset(f_Urad).sfcWind.sel(time=slice('0206-01', '0235-12')) * (4.87 / log(67.8 * 10 - 5.42))
# f_Ubgc = r"E:\CMIP6\ACCESS\sfcWind_Amon_ACCESS-ESM1-5_1pctCO2-bgc_r1i1p1f1_gn_010101-025012.nc"
# Ubgcq30 = xr.open_dataset(f_Ubgc).sfcWind.sel(time=slice('0102-01', '0131-12')) * (4.87 / log(67.8 * 10 - 5.42))
# Ubgch30 = xr.open_dataset(f_Ubgc).sfcWind.sel(time=slice('0206-01', '0235-12')) * (4.87 / log(67.8 * 10 - 5.42))

# UCTL = UCTLh30.mean(dim='time') - UCTLq30.mean(dim='time')
# Urad = Uradh30.mean(dim='time') - Uradq30.mean(dim='time')
# Ubgc = Ubgch30.mean(dim='time') - Ubgcq30.mean(dim='time')

# xianghu = UCTL - Urad - Ubgc



data = adjust_longitude(xianghu)
rdata = bilinear_interpolation(data)
rdata.to_netcdf(r"E:\CMIP6\UKESM\RN_int.nc")

