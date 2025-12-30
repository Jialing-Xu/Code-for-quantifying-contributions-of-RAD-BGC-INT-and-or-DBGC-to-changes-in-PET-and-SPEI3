import numpy as np
import xarray as xr
from climate_indices import indices
from climate_indices import compute  # 计算SPEI的包

f_prrad = r"E:\mingan_results\ACCESS\1P_no_rad.nc"
prerad = xr.open_dataset(f_prrad).pr
f_prbgc = r"E:\mingan_results\ACCESS\1P_no_bgc.nc"
prebgc = xr.open_dataset(f_prbgc).pr
f_pr3 = r"E:\mingan_results\ACCESS\1P_no_int.nc"
preint = xr.open_dataset(f_pr3).pr

f_pr = r"E:\CMIP6\ACCESS\pr_1pctCO2_1850-1999.nc"
pre1 = xr.open_dataset(f_pr).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24
pre = xr.open_dataset(f_pr).pr.sel(time=slice('1850-01', '1984-12')) *60*60*24
pre_54 = xr.open_dataset(f_pr).pr.sel(time=slice('1954-01', '1954-12')) *60*60*24

f_tas1 = "E:\\PET_results\\ACCESS\\1PET_tas_1.nc"
pet_tas1 = xr.open_dataset(f_tas1).PET
f_tas2 = "E:\\PET_results\\ACCESS\\1PET_tas_2.nc"
pet_tas2 = xr.open_dataset(f_tas2).PET
f_tas3 = "E:\\PET_results\\ACCESS\\1PET_tas_3.nc"
pet_tas3 = xr.open_dataset(f_tas3).PET
f_rn1 = "E:\\PET_results\\ACCESS\\1PET_rn_1.nc"
pet_rn1 = xr.open_dataset(f_rn1).PET
f_rn2 = "E:\\PET_results\\ACCESS\\1PET_rn_2.nc"
pet_rn2 = xr.open_dataset(f_rn2).PET
f_rn3 = "E:\\PET_results\\ACCESS\\1PET_rn_3.nc"
pet_rn3 = xr.open_dataset(f_rn3).PET
f_rh1 = "E:\\PET_results\\ACCESS\\1PET_rh_1.nc"
pet_rh1 = xr.open_dataset(f_rh1).PET
f_rh2 = "E:\\PET_results\\ACCESS\\1PET_rh_2.nc"
pet_rh2 = xr.open_dataset(f_rh2).PET
f_rh3 = "E:\\PET_results\\ACCESS\\1PET_rh_3.nc"
pet_rh3 = xr.open_dataset(f_rh3).PET
f_u21 = "E:\\PET_results\\ACCESS\\1PET_u2_1.nc"
pet_u21 = xr.open_dataset(f_u21).PET
f_u22 = "E:\\PET_results\\ACCESS\\1PET_u2_2.nc"
pet_u22 = xr.open_dataset(f_u22).PET
f_u23 = "E:\\PET_results\\ACCESS\\1PET_u2_3.nc"
pet_u23 = xr.open_dataset(f_u23).PET
f_co2 = "E:\\PET_results\\ACCESS\\PET_co2.nc"
pet_co2 = xr.open_dataset(f_co2).PET

f_pet = r"E:\PET_results\ACCESS\CTL_PET.nc"
pet1 = xr.open_dataset(f_pet).PET.sel(time=slice('1955-01', '1984-12'))
pet = xr.open_dataset(f_pet).PET.sel(time=slice('1850-01', '1984-12'))
pet_54 = xr.open_dataset(f_pet).PET.sel(time=slice('1954-01', '1954-12'))

# 拼接数据
combined_pre = xr.concat([pre, pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, pre1,
                          pre_54, prerad,
                          pre_54, prebgc,
                          pre_54, preint,
                          pre_54, pre1
                          ], dim='time', coords='minimal', compat='override')
combined_pet = xr.concat([pet, pet_54, pet_tas1,
                          pet_54, pet_tas2,
                          pet_54, pet_tas3,
                          pet_54, pet_rn1,
                          pet_54, pet_rn2,
                          pet_54, pet_rn3,
                          pet_54, pet_rh1,
                          pet_54, pet_rh2,
                          pet_54, pet_rh3,
                          pet_54, pet_u21,
                          pet_54, pet_u22,
                          pet_54, pet_u23,
                          pet_54, pet1,
                          pet_54, pet1,
                          pet_54, pet1,
                          pet_54, pet_co2
                          ], dim='time', coords='minimal', compat='override')

# 计算每个月的天数
month_lengths = combined_pre.time.dt.days_in_month

# 计算每个月总的数据
pre_total = combined_pre * month_lengths

# 数据读取转换为np.array
pre_data = np.asarray(pre_total)
pet_data = np.asarray(combined_pet)
pre_data[pre_data < 0] = 0

# 提取时间、经度、纬度信息
time = combined_pre['time']
longitude = combined_pre['lon']
latitude = combined_pre['lat']

# 计算SPEI
spei = np.zeros((7572, 145, 192), dtype=float, order='C')*np.nan

for lat_index in range(len(latitude)):
    for lon_index in range(len(longitude)):
        pet_data_grid = pet_data[:, lat_index, lon_index]
        pre_data_grid = pre_data[:, lat_index, lon_index]

        spei[:, lat_index, lon_index] = indices.spei(precips_mm=pre_data_grid,
                                                     pet_mm=pet_data_grid,
                                                     scale=3,  # 3个月尺度
                                                     distribution=indices.Distribution,
                                                     periodicity=compute.Periodicity.monthly,
                                                     data_start_year=1851,
                                                     calibration_year_initial=1851,
                                                     calibration_year_final=1984,
                                                     )


data = xr.Dataset(
    {
        "spei": (["time", "lat", "lon"], spei),
    },
    coords={
        "time": combined_pre.coords["time"],
        "lat": combined_pre.coords["lat"],
        "lon": combined_pre.coords["lon"],
    },
)

output_file = r"E:\SPEI_results\ACCESS\1all_SPEI3.nc"
data.to_netcdf(output_file)

print(spei)

