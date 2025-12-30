import xarray as xr
from math import pi, sin, cos, tan, log

# ##### CO2前三十年
# f_co2CTL = r"E:\CMIP6\UKESM\co2s_1%.nc"
# co2CTLq30 = xr.open_dataset(f_co2CTL).co2s.sel(time=slice('1851-01', '1880-12'))

# co2CTL = co2CTLq30.groupby('time.month').mean('time')

# months = co2CTLq30['time.month'].values
# co2q30 = co2CTL.sel(month=months).values

# data1 = xr.Dataset(
#     {"co2s": (["time", "lat", "lon"], co2q30)},
#     coords={
#         "time": co2CTLq30.coords["time"],
#         "lat": co2CTLq30.coords["lat"],
#         "lon": co2CTLq30.coords["lon"],
#     }
# )

# data1.to_netcdf(r"E:\mingan_results\UKESM\co2_q30.nc")




# ##### 降水变化
# f_PCTL = r"E:\CMIP6\UKESM\pr_1pctCO2_1850-1999.nc"
# PCTLq30 = xr.open_dataset(f_PCTL).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
# PCTLh30 = xr.open_dataset(f_PCTL).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24
# f_Prad = r"E:\CMIP6\UKESM\pr_1pctCO2-rad_1850-1999.nc"
# Pradq30 = xr.open_dataset(f_Prad).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
# Pradh30 = xr.open_dataset(f_Prad).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24
# f_Pbgc = r"E:\CMIP6\UKESM\pr_1pctCO2-bgc_1850-1999.nc"
# Pbgcq30 = xr.open_dataset(f_Pbgc).pr.sel(time=slice('1851-01', '1880-12')) *60*60*24
# Pbgch30 = xr.open_dataset(f_Pbgc).pr.sel(time=slice('1955-01', '1984-12')) *60*60*24

# months = PCTLh30['time.month'].values

# PCTL = PCTLh30 - PCTLq30.groupby('time.month').mean('time').sel(month=months).values
# Prad = Pradh30 - Pradq30.groupby('time.month').mean('time').sel(month=months).values
# Pbgc = Pbgch30 - Pbgcq30.groupby('time.month').mean('time').sel(month=months).values
# Pint = PCTL - Prad - Pbgc

# no_rad = PCTLh30 - Prad
# no_bgc = PCTLh30 - Pbgc
# no_int = PCTLh30 - Pint

# no_rad.to_netcdf(r"E:\mingan_results\UKESM\1P_no_rad.nc")
# no_bgc.to_netcdf(r"E:\mingan_results\UKESM\1P_no_bgc.nc")
# no_int.to_netcdf(r"E:\mingan_results\UKESM\1P_no_int.nc")



# ##### 温度变化
# f_TmaxCTL = r"E:\CMIP6\UKESM\tasmax_1pctCO2_1850-1999.nc"
# TmaxCTLq30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('1851-01', '1880-12')) - 273.16
# TmaxCTLh30 = xr.open_dataset(f_TmaxCTL).tasmax.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tmaxrad = r"E:\CMIP6\UKESM\tasmax_1pctCO2-rad_1850-1999.nc"
# Tmaxradq30 = xr.open_dataset(f_Tmaxrad).tasmax.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tmaxradh30 = xr.open_dataset(f_Tmaxrad).tasmax.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tmaxbgc = r"E:\CMIP6\UKESM\tasmax_1pctCO2-bgc_1850-1999.nc"
# Tmaxbgcq30 = xr.open_dataset(f_Tmaxbgc).tasmax.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tmaxbgch30 = xr.open_dataset(f_Tmaxbgc).tasmax.sel(time=slice('1955-01', '1984-12')) - 273.16

# months = TmaxCTLh30['time.month'].values

# TmaxCTL = TmaxCTLh30 - TmaxCTLq30.groupby('time.month').mean('time').sel(month=months).values
# Tmaxrad = Tmaxradh30 - Tmaxradq30.groupby('time.month').mean('time').sel(month=months).values
# Tmaxbgc = Tmaxbgch30 - Tmaxbgcq30.groupby('time.month').mean('time').sel(month=months).values
# Tmaxint = TmaxCTL - Tmaxrad - Tmaxbgc

# no_rad = TmaxCTLh30 - Tmaxrad
# no_bgc = TmaxCTLh30 - Tmaxbgc
# no_int = TmaxCTLh30 - Tmaxint

# no_rad.to_netcdf(r"E:\mingan_results\UKESM\1Tmax_no_rad.nc")
# no_bgc.to_netcdf(r"E:\mingan_results\UKESM\1Tmax_no_bgc.nc")
# no_int.to_netcdf(r"E:\mingan_results\UKESM\1Tmax_no_int.nc")


# f_TminCTL = r"E:\CMIP6\UKESM\tasmin_1pctCO2_1850-1999.nc"
# TminCTLq30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('1851-01', '1880-12')) - 273.16
# TminCTLh30 = xr.open_dataset(f_TminCTL).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tminrad = r"E:\CMIP6\UKESM\tasmin_1pctCO2-rad_1850-1999.nc"
# Tminradq30 = xr.open_dataset(f_Tminrad).tasmin.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tminradh30 = xr.open_dataset(f_Tminrad).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16
# f_Tminbgc = r"E:\CMIP6\UKESM\tasmin_1pctCO2-bgc_1850-1999.nc"
# Tminbgcq30 = xr.open_dataset(f_Tminbgc).tasmin.sel(time=slice('1851-01', '1880-12')) - 273.16
# Tminbgch30 = xr.open_dataset(f_Tminbgc).tasmin.sel(time=slice('1955-01', '1984-12')) - 273.16

# months = TminCTLh30['time.month'].values

# TminCTL = TminCTLh30 - TminCTLq30.groupby('time.month').mean('time').sel(month=months).values
# Tminrad = Tminradh30 - Tminradq30.groupby('time.month').mean('time').sel(month=months).values
# Tminbgc = Tminbgch30 - Tminbgcq30.groupby('time.month').mean('time').sel(month=months).values
# Tminint = TminCTL - Tminrad - Tminbgc

# no_rad = TminCTLh30 - Tminrad
# no_bgc = TminCTLh30 - Tminbgc
# no_int = TminCTLh30 - Tminint

# no_rad.to_netcdf(r"E:\mingan_results\UKESM\1Tmin_no_rad.nc")
# no_bgc.to_netcdf(r"E:\mingan_results\UKESM\1Tmin_no_bgc.nc")
# no_int.to_netcdf(r"E:\mingan_results\UKESM\1Tmin_no_int.nc")



# ##### 净辐射变化
# f_rssCTL = r"E:\CMIP6\UKESM\rss_1pctCO2_1850-1999.nc"
# RnsCTLq30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1851-01', '1880-12'))
# RnsCTLh30 = xr.open_dataset(f_rssCTL).rss.sel(time=slice('1955-01', '1984-12'))
# f_rlsCTL = r"E:\CMIP6\UKESM\rls_1pctCO2_1850-1999.nc"
# RnlCTLq30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1851-01', '1880-12'))
# RnlCTLh30 = xr.open_dataset(f_rlsCTL).rls.sel(time=slice('1955-01', '1984-12'))

# RNCTLq30 = (RnsCTLq30 + RnlCTLq30) * 0.0864
# RNCTLh30 = (RnsCTLh30 + RnlCTLh30) * 0.0864

# f_rssrad = r"E:\CMIP6\UKESM\rss_1pctCO2-rad_1850-1999.nc"
# Rnsradq30 = xr.open_dataset(f_rssrad).rss.sel(time=slice('1851-01', '1880-12'))
# Rnsradh30 = xr.open_dataset(f_rssrad).rss.sel(time=slice('1955-01', '1984-12'))
# f_rlsrad = r"E:\CMIP6\UKESM\rls_1pctCO2-rad_1850-1999.nc"
# Rnlradq30 = xr.open_dataset(f_rlsrad).rls.sel(time=slice('1851-01', '1880-12'))
# Rnlradh30 = xr.open_dataset(f_rlsrad).rls.sel(time=slice('1955-01', '1984-12'))

# RNradq30 = (Rnsradq30 + Rnlradq30) * 0.0864
# RNradh30 = (Rnsradh30 + Rnlradh30) * 0.0864

# f_rssbgc = r"E:\CMIP6\UKESM\rss_1pctCO2-bgc_1850-1999.nc"
# Rnsbgcq30 = xr.open_dataset(f_rssbgc).rss.sel(time=slice('1851-01', '1880-12'))
# Rnsbgch30 = xr.open_dataset(f_rssbgc).rss.sel(time=slice('1955-01', '1984-12'))
# f_rlsbgc = r"E:\CMIP6\UKESM\rls_1pctCO2-bgc_1850-1999.nc"
# Rnlbgcq30 = xr.open_dataset(f_rlsbgc).rls.sel(time=slice('1851-01', '1880-12'))
# Rnlbgch30 = xr.open_dataset(f_rlsbgc).rls.sel(time=slice('1955-01', '1984-12'))

# RNbgcq30 = (Rnsbgcq30 + Rnlbgcq30) * 0.0864
# RNbgch30 = (Rnsbgch30 + Rnlbgch30) * 0.0864

# months = RNCTLh30['time.month'].values

# RNCTL = RNCTLh30 - RNCTLq30.groupby('time.month').mean('time').sel(month=months).values
# RNrad = RNradh30 - RNradq30.groupby('time.month').mean('time').sel(month=months).values
# RNbgc = RNbgch30 - RNbgcq30.groupby('time.month').mean('time').sel(month=months).values
# RNint = RNCTL - RNrad - RNbgc

# no_rad = (RNCTLh30 - RNrad).where(lambda x: x >= 0, 0)
# no_bgc = (RNCTLh30 - RNbgc).where(lambda x: x >= 0, 0)
# no_int = (RNCTLh30 - RNint).where(lambda x: x >= 0, 0)

# no_rad.to_netcdf(r"E:\mingan_results\UKESM\1RN_no_rad.nc")
# no_bgc.to_netcdf(r"E:\mingan_results\UKESM\1RN_no_bgc.nc")
# no_int.to_netcdf(r"E:\mingan_results\UKESM\1RN_no_int.nc")



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

# months = RHCTLh30['time.month'].values

# RHCTL = RHCTLh30 - RHCTLq30.groupby('time.month').mean('time').sel(month=months).values
# RHrad = RHradh30 - RHradq30.groupby('time.month').mean('time').sel(month=months).values
# RHbgc = RHbgch30 - RHbgcq30.groupby('time.month').mean('time').sel(month=months).values
# RHint = RHCTL - RHrad - RHbgc

# no_rad = RHCTLh30 - RHrad
# no_bgc = RHCTLh30 - RHbgc
# no_int = RHCTLh30 - RHint

# no_rad.to_netcdf(r"E:\mingan_results\UKESM\1RH_no_rad.nc")
# no_bgc.to_netcdf(r"E:\mingan_results\UKESM\1RH_no_bgc.nc")
# no_int.to_netcdf(r"E:\mingan_results\UKESM\1RH_no_int.nc")



##### 风速变化
f_UCTL = r"E:\CMIP6\UKESM\sfcWind_1pctCO2_1850-1999.nc"
UCTLq30 = xr.open_dataset(f_UCTL).sfcWind.sel(time=slice('1851-01', '1880-12')) * (4.87 / log(67.8 * 10 - 5.42))
UCTLh30 = xr.open_dataset(f_UCTL).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))
f_Urad = r"E:\CMIP6\UKESM\sfcWind_1pctCO2-rad_1850-1999.nc"
Uradq30 = xr.open_dataset(f_Urad).sfcWind.sel(time=slice('1851-01', '1880-12')) * (4.87 / log(67.8 * 10 - 5.42))
Uradh30 = xr.open_dataset(f_Urad).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))
f_Ubgc = r"E:\CMIP6\UKESM\sfcWind_1pctCO2-bgc_1850-1999.nc"
Ubgcq30 = xr.open_dataset(f_Ubgc).sfcWind.sel(time=slice('1851-01', '1880-12')) * (4.87 / log(67.8 * 10 - 5.42))
Ubgch30 = xr.open_dataset(f_Ubgc).sfcWind.sel(time=slice('1955-01', '1984-12')) * (4.87 / log(67.8 * 10 - 5.42))

months = UCTLh30['time.month'].values

UCTL = UCTLh30 - UCTLq30.groupby('time.month').mean('time').sel(month=months).values
Urad = Uradh30 - Uradq30.groupby('time.month').mean('time').sel(month=months).values
Ubgc = Ubgch30 - Ubgcq30.groupby('time.month').mean('time').sel(month=months).values
Uint = UCTL - Urad - Ubgc

no_rad = UCTLh30 - Urad
no_bgc = UCTLh30 - Ubgc
no_int = UCTLh30 - Uint

no_rad.to_netcdf(r"E:\mingan_results\UKESM\1U_no_rad.nc")
no_bgc.to_netcdf(r"E:\mingan_results\UKESM\1U_no_bgc.nc")
no_int.to_netcdf(r"E:\mingan_results\UKESM\1U_no_int.nc")



