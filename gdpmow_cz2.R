
##
## Nowcast Petrova hodnota
##
#------------------------------------------------------------------------------------------------------------------------------------------------------#
# load libraries----
#------------------------------------------------------------------------------------------------------------------------------------------------------#
rm(list=ls())

options(prompt = "R> ")




source("Eurostat_functions.R")
source("other_functions.R")
source("funkce_rady.R")
library("xts")
require("dplyr")
require("reshape2")
require("forecast")
require("ggplot2")
require("tseries")
require("rsdmx")
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/NOWCAST/NowcastModSelect.R")
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/NOWCAST/export_nowcast_results.R")
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/TIMESERIES/xts_ts_conversions.R")
require("dygraphs")
library("eurostat")
library("rvest")
library("knitr")
library("highr")
require(data.table)
require("plyr")

## Estimation & measuring forecasting performance tools
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/rfvar_m.R") # Reduced form BVAR estimation
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/postdraw_m.R") # draws from BVAR posterior distribution
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/rwwish.R") # draws from Wishart distribution
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/rballunif.R") # draws from unit sphere
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/fcast.R") # forecasting
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/mgnldnsty_m.R") # Reduced form BVAR estimation - with Minnesota prior
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/matrictint.R") # Scale factor for a matrix t distribution, like the posterior from a VAR
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/varprior.R") # Minnesota prior via dummy observations
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/stability_test.R") # test of dynamic stability of VAR
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/bvar_performance.R") # test of performance of Minnesota prior BVAR

## Structural analysis tools
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/impulsdtrf.R") # impulse response functions
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/Uhlig_reject_m.R") # IRF via sign restrictions - Uhlig reject method
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/Uhlig_accept.R") # IRF via sign restrictions - Uhlig acceptation algorithm
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/BVAR_Tools/yield_curve.R")#yield curve factors
## Libraries
library("xts")
library("lubridate")
require("forecast")
require("MacrobondAPI")

## Other functions
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/FAPPLY/FAPPLY.R") # data processing tools
source("C:/Users/xburj/OneDrive/Dokumenty/R_projects/functions/csob_colors.R") # CSOB colors for graphs
minf <- function(x) mean(x, na.rm = T)
export_nowcast_results <- function(model, nowcast, indicator = NULL, plot_title = NULL, plot_subtitle = NULL, period = NULL) {
  
  if(class(model) != "bsts") stop("'model' has to be of class 'bsts'!")
  if(class(nowcast) != "bsts.prediction") stop("'nowcast' has to be of class 'bsts.prediction'!")
  if(is.null(indicator)) stop("'indicator' is missing!")
  if(is.null(plot_title)) stop("'plot_title' is missing!")
  if(is.null(plot_subtitle)) stop("'plot_subtitle' is missing!")
  
  m <- model[c("coefficients", "state.specification", "original.series", "predictors", "niter")]
  n <- nowcast$distribution
  dt <- Sys.Date()
  
  res <- list(date = dt, plot_title = plot_title, plot_subtitle = plot_subtitle, model = m, nowcast = n, period = period)
  
  dir_path <- paste0("u:/Dokumenty/Rdata/GDPnow_emu/", indicator, "/", period)
  full_path <- paste0(dir_path, "/", gsub("-", "", dt), ".rds")
  
  if(!dir.exists(dir_path)) dir.create(dir_path, recursive = T)
  
  saveRDS(res, full_path)
  message("--- Nowcast data exported to: ", full_path, " ---", sep = "")
  
}


#------------------------------------------------------------------------------------------------------------------------------------------------------#
# initial assumptions----
#------------------------------------------------------------------------------------------------------------------------------------------------------#

st.time <- 2000
en.time <- 2017

logdif <- c("indu_emu",   "manu_emu",   "coreretail_emu",    "retail_de",  "retail_fr",  "retail_it",  "manu_de",    "manu_fr",    "manu_it",    "const_emu",  "retail_emu", "trade_emu","unem_emu","hourwork_de")
depend  <- "gdp_emu"
regresor <- c("indu_emu",   "manu_emu",   "coreretail_emu", "const_emu",  "retail_emu", "trade_emu","unem_emu","hourwork_de", "ifo_de","pastprod_indu", "confidence_indu","employ_indu", "demand_serv")
aaa <- "2022-12-01" # the last period to be forecasted from
bbb <- "2023-01-01" # period to be forecasted (vzdy beyprostredne po)
regresor_f  <- c("manu_emu","retail_emu", "const_emu") # final selection of the regressoprs for the model
st.date_graph  <- c(2012,1)
l <- 0.01 #jake promenne pustit do modelu

#------------------------------------------------------------------------------------------------------------------------------------------------------#

# Download data-----

#------------------------------------------------------------------------------------------------------------------------------------------------------#
#Data from eurostat
{
  toc <- get_eurostat_toc()
  kable(head(toc))
  
  gdp <- get_eurostat(id = "namq_10_gdp", time_format = "num")
  gdp<- as.data.table(gdp)
  gdp <- gdp[geo =="CZ" & unit %in% c("CLV_PCH_PRE","CLV05_MNAC")&s_adj=="SCA"&na_item=="B1GQ",
             .(unit,time,values)]
  gdp <- dcast.data.table(data = gdp,formula = time~unit,fun.aggregate = mean)
  gdp_xts <- xts(x =gdp[,c("CLV05_MNAC","CLV_PCH_PRE")] ,order.by = seq(from=as.Date(paste0(first(gdp$time),"-03-01")),by = "3 months",length.out = nrow(gdp)))
  #gdp_select <- gdp_xts[,"CLV_PCH_PRE"]
  gdp_select <- gdp_xts
  colnames(gdp_select) <- c( "gdp_vol","gdp_emu")
  
    manu <- get_eurostat(id = "sts_inpr_m", time_format = "num")
    manu<- as.data.table(manu)
  manu <- manu[indic_bt=="PROD"&s_adj=="SCA"&unit=="I15"&geo=="CZ",
               .(nace_r2,time,values)]
  manu <- dcast.data.table(data = manu,formula = time~nace_r2,fun.aggregate = mean)
  manu_xts <- xts(x = manu[,-1],order.by = seq(from=as.Date(paste0(first(manu$time),"-01-01")),by = "1 months",length.out = nrow(manu)))
  manu_xts <- (diff.xts(log(apply.quarterly(manu_xts,FUN = mean)),lag = 1))*100
  
  
  retail <- get_eurostat(id = "sts_trtu_m", time_format = "num")
  retail <- as.data.table(retail)
  retail <- retail[indic_bt=="TOVV"&s_adj=="SCA"&unit=="I15"&geo=="CZ",
                   .(nace_r2,time,values)]
  retail <- dcast.data.table(data = retail,formula = time~nace_r2,fun.aggregate = mean)
  retail_xts <- xts(x = retail[,-1],order.by = seq(from=as.Date(paste0(first(retail$time),"-01-01")),by = "1 months",length.out = nrow(retail)))
  retail_xts <- (diff.xts(log(apply.quarterly(retail_xts,FUN = mean)),lag = 1))*100
  
  service <- get_eurostat(id = "sts_setu_m", time_format = "num")
  service <- as.data.table(service)
  unique(service$unit)
  service <- service[indic_bt=="TOVT"&s_adj=="SCA"&unit=="I15"&geo=="DE",
                     .(nace_r2,time,values)]
  service <- dcast.data.table(data = service,formula = time~nace_r2,fun.aggregate = mean)
  
  service_xts <- xts(x = service[,-1],order.by = seq(from=as.Date(paste0(first(service$time),"-01-01")),by = "1 months",length.out = nrow(service)))
  service_xts <- (diff.xts(log(apply.quarterly(service_xts,FUN = mean)),lag = 1))*100
  
  
  const <- get_eurostat(id = "sts_copr_m", time_format = "num")
  const <- as.data.table(const)
  const <- const[indic_bt=="PROD"& s_adj=="SCA"&unit=="I15"&geo=="CZ",
                 .(nace_r2,time,values)]
  const <- dcast.data.table(data = const,formula = time~nace_r2,fun.aggregate = mean)
  const_xts <- xts(x = const[,-1],order.by = seq(from=as.Date(paste0(first(const$time),"-01-01")),by = "1 months",length.out = nrow(const)))
  const_xts <- (diff.xts(log(apply.quarterly(const_xts,FUN = mean)),lag = 1))*100
  
  unem <- get_eurostat(id = "une_rt_m", time_format = "num")
  unem <- as.data.table(unem)
  unem <- unem[s_adj=="SA"&age=="TOTAL"&sex=="T"&geo=="CZ"&unit=="THS_PER",
               .(time,values)]
  setkey(x = unem,time)
  unem_xts <- xts(x = unem[,-1],order.by =seq(from=as.Date(paste0(substr(first(unem$time),start = 0,stop = 4),"-04-01")),by = "1 months",length.out = nrow(unem)))
  unem_xts <- (diff.xts(log(apply.quarterly(unem_xts,FUN = mean)),lag = 1))*100
  colnames(unem_xts) <- "unem"
  
  #int_trade <- get_eurostat(id = "ext_st_28msbec", time_format = "num")
  #int_trade <- as.data.table(int_trade)
  #int_trade <- int_trade[indic_et=="TRD_VAL_SCA"&partner=="WORLD"&bclas_bec=="TOTAL"&geo=="CZ"&stk_flow%in%c("EXP","IMP"),
                         #.(time,values,stk_flow)]
  #setkey(x = int_trade,time)
  #int_trade <- dcast.data.table(data = int_trade,formula = time~stk_flow,fun.aggregate = mean,value.var = "values")
  
  #int_trade_xts <- xts(x = int_trade[,-1],order.by =seq(from=as.Date(paste0(first(int_trade$time),"-01-01")),by = "1 months",length.out = nrow(int_trade)))
  #int_trade_xts <- (diff.xts(log(apply.quarterly(int_trade_xts,FUN = mean)),lag = 1))*100
  
  #colnames(int_trade_xts) <- c("exp","imp")
  
  

  #dat_q <- readRDS("dat_q.rds")
  #dat_m <- readRDS("dat_m.rds")
}

#Data from Macrobond
{
  #Download data
  #tic <- list(US = list(GDP = "usnaac0169", U_RATE = "uslama1849", INFL = "uspric2156", PCE_core = "uspric0006", BRENT = "uscaes0302", BAA_Y = "usrate0217", US_10Y = "us10ygov", VIX = "vix", EURUSD = "usdecbfix", FED_MB = "usmost0111", FED_FUNDS = "usrate0190", UScreditBIS = "bis_quspamxdca", EURcreditBIS= "bis_qxmpamxdca", EURloans2 = "eubank0430"),
  #            EU = list(GDP = "eunaac2903", U_RATE = "eulama0019", INFL_CORE = "eupric0151", IT_10Y = "it10ygov", DE_10Y = "de10ygov", ECB_ASSETS = "eubank0133", EONIA = "deonibr"))
  
  tic <- list(US = list(GDPrealChained = "usnaac0169", GDPcurrentChained = "usnaac0057",
                        GDPindex = "usnaac0277", GDPrealpotential = "usfcst0031", PersonalSpendingSA = "usnaac0591", RealRetailSales = "ustrad1000", Payrolls_All = "uslama1060", Payrolls_private = "uslama1061", Initial_claims = "uslama3295", Continuing_claims = "uslama3350",  Covered_employment = "uslama4757", ISM_man = "ussurv1055", ISM_services = "ussurv1044", CapUtil_Industry = "usprod1071", PCEcoreIndex = "uspric0006", PCEIndex = "uspric0001", CPI = "uspric2156",  coreCPI = "uspric2373", uRate = "uslama1849", uRate_Manufacturing = "bls_lnu04032232m", FedFunds = "usrate0190", TreasuryRate3M = "us03mgov", UST10Y = "us10Ygov", M2money = "usmost0002", FedMonetaryBase = "usmost0111", FedM1SA = "usmost0001", FedM2SA= "usmost0002", Bonds_at_Fed = "usbank0144", UScomloans = "usbank0685", USrestateloans = "usbank0686", ConsumerCreditSA= "usbank0103", ConsumerCreditWeeklySA = "usbank0498", HourlyEarnSA = "uslama0218", AtlFed_wage_tracker = "usinea0234", CA_balance = "usbopa0724", Labor_force = "bls_lns11000000",
                        Output_per_hourSA = "uslama0202", Output_per_hourSAAR = "uslama0201", Output_per_hourYY = "uslama0200", STLFSI = "ussurv0100", CNFCI = "ussurv0386"),
              EA = list(GDPreal = "eunaac0400", GDPnom = "eunaac0001", GDPrealpotential_EC = "ameco_ea19_1_0_0_0_ovgdp",GDPnomSA = "eunaac0019",
                        GDPindex = "buba_mb_70070", OutputGap = "eufcst0024",
                        CapUtil_Ind_quarter = "euprod0496", IndustryProd_ECB = "ecb_00213147", PMI_Man = "ih:mb:com:pmi_emu_manufacturing", PMI_Com = "ih:mb:com:pmi_emu_composite",
                        HICPnsa= "eupric0001", HICPsa= "ecb_00678029", HICPcore = "i05totxnrgfoodeadx", HICPcoreSA = "ecb_00678040", German_goods_CPI = "depric2453", uRate = "eulama0019", ULC = "ecb_00171346", Wages_negotiatedECB = "euinea0071", FranceRate3M = "fr3mgov", EONIA = "deonibr", M1moneySA = "eumost0062", M2moneySA = "eumost0088", M3moneySA = "eumost0110", ECBassets = "eubank0133", EURloans1 = "eubank0050", CompensationSA = "eunaac0034", EmploymentNSA = "eulama3457", UnemploymentSA = "eulama0001",ifo="desurv1001", ECBsurvey_infl_exp1Y = "eufcst0007", ECBsurvey_infl_exp2Y = "eufcst0008", ECBsurvey_infl_expLR = "eufcst0009", ECB_system_Stress_Index = "ecb_cissdu2z0z4fecss_ciidx"),
              FX = list(EURUSD = "usdecbfix", USDindexReal = "usexri0005", USDindexNominal = "usexri0008",
                        CZK_REER_HICP = "czecexri0002", CZK_REER_PPI = "czexri0003", CZK_REER_ULC = "czexri0007", CZK_REER_Defl = "czexri0011", EUR_REERecb38  = "ecb_00170141",
                        EUR_EERecb38 = "ecb_00137743", CNY_broad = "bis_dnbcn", EURPLN = "plnecbfix", USDPLN = "pln", USDTRY = "try"),
              CM = list(WorldComIndexHWWI = "wocind0025", CRBindex = "wocind0250", Brent = "uscaes0302", WTI = "wocaes0076", Wheat = "ebm_c1_cl"),
              EX = list(BaaSpread = "usrate0217", VIX = "vix", SP500 = "sp500_500", SP500real = "sp500_realpr", DAX = "deeqin0032", EuroStoXX_Vol3M = "vstoxx3m", IT_10Y = "it10ygov", DE_10Y = "de10ygov", USLaborForceSA = "bls_lns11000000", USLaborForceSA16over = "uslama1698", USpopulation15_74SA = "oecd_mei_usa_lfwa74tt_stsa_m"),
              CZ = list(czcorehicpi = "i15totxnrgfoodczdx",PMI_CZ="ih:mb:com:pmi_man_cz",czgdp="cznaac0246", GDPrealChngQoQ = "cznaac0218", czcorploan="czbank0134",
                        czhicp="i15cp00czdx",czvacan="czlama0111",czcpi="czpric0001",eurczk="eurcnbczkfixaver",czgdpind="clvi10nsab1gqcz1g",empl_rate="czlama0002",
                        czunem="czlama0051", uRate_sa = "czlama0060",ind_prod="czprod0041",new_car="ecb_stsmczycregpc00003abs",retail="cztrad0562",trade_balance="cztrad0005",
                        terms_trad="czpric0177", bus_surv="czsurv0006",empl_m="czlama0066",cons_credit="czbank0178",hous_lend="czbank0180",
                        machinery="czprod0059",car_prod="czprod0060",steel_prod="czsteelprod",ec_recentprod="czecfin0003",
                        euro_induconf="bs_ici_balnsaczib",ec_ecsent="czecfin0001",ec_bussent="czecfin0002",czso_ecsent="czsurv0008",
                        transport="cztrad0390",office_sup="cztrad0510",car_balance="cztrad0193", machinery_exp="cztrad0122",cz_warehouse="cztrad0406",
                        cz_freight="cztrad0398",cz_advertise="cztrad0484",emu_constr="bs_cci_balnsaczib",emu_consum="bs_csmci_balnsaczib",
                        emu_indu="bs_ici_balnsaczib",ec_constrempl="czecfin0100",ec_machinorder="indu_cz_28_2_bs_m",ec_autoorder="indu_cz_29_2_bs_m",
                        ec_invorder="indu_cz_intm_2_bs_m",ec_induempl="czecfin0099", house_rent_oecd="oecd_eo_hpi_00463236",
                        house_price_oecd="oecd_eo_hpi_00463387",house_princome_oecd="oecd_eo_hpi_00463265",house_prrent_oecd="oecd_eo_hpi_00463318", EU_funds_transfers = "mioeurkas1s1creeuicz6q",electricity_cons="ih:mb:com:elec_load_wh_cz"),
              PL = list(GDPrealChained = "plnaac0197", GDPreal_yy_GUS = "plnaac0017", GDPreal_qq_GUS = "plnaac0065", GDPcurrentChained = "cpmnacscab1gqpl1g", Productivity_OECD = "oecd_mei_00423591",
                        FixedInvestments = "clv10mnacscap51gpl1g", ConsumptionSA_GUS = "plnaac0190", ConsumptionPerCap = "clvpchsmhabnsap31s14s15plmq", NumWorkDays = "pllama5000", Vacancies = "pllama0145",
                        VacanciesUnfilled = "pllama0144", EmploymentTotal = "pllama0004", UnemploymentTotalReg = "pllama0132", uRateReg = "pllama0140", uRate_SA = "sattotalpcactplq9",
                        EmploymentEnterprises = "pllama0050", EmploymentEnterprisesAVG = "pllama0167", WagesPrivate = "plinea0001", WageReal_yy = "plinea0051",
                        RetailTradeRealyy = "pltrad0103", RetailTradeRealSA = "pltrad0151", CarRegistrationSA_ECB = "ecb_stsmplycregpc00003abs", RetailConfidenceSA_EC = "plecfin0070",
                        IndustProdSA = " plprod0063", IndustProd_yy = "plprod0065", ConstructionSA = "plprod0180", Construction_yy = "plprod0181", CPI= "plpric0005", HICP = "plpric0573",
                        coreCPIyy = "plpric0303", CPIyyAdminExcl = "plpric0301", CPIyyFood = "plpric0023", CPIyyBread = "plpric0493", Wheat_price = "plpric0171", PPIyy = "plpric0065",
                        PPI_clothing_yy = "plpric0073", PPI_autos_yy = "plpric0088", PPI_furniture_yy = "plpric0090", GDPexternalTradePLN = "cpmnacnsab11pl1g", ExportGoods_EUR_NBP = "plbopa0003",
                        ImportGoods_EUR_NBP = "plbopa0004", CapitalAccountBalEUR = "plbopa0457", CapitalAccountCreditPLN  = "plbopa0199", CapitalAccountCredit_EUR = "plbopa0015",
                        TradeBalanceGUS = "pltrad0102", EU_funds_transfers = "mioeurkas1s1creeuipl6q", REER_ECB = "ecb_00170225", REER_IMF = "imfm964ereer_ix", M3 = "plmost0004", FX_ReservesEUR = "plfofi0040"),
              DE = list(PMI_Man="ih:mb:com:pmi_man_de" , IFO_MANU="desurv0022",IFO_CONST="desurv0023",IFO_TRADE="deifo_ghanges0_kld_bds",IFO_SERV="desurv0341",IFO_DE="desurv01046", IFO_Man_expct = "desurv1011",
                        IFO_Man_Total = "desurv1007", vacancies_DE="delama0016",manufact="deprod1414",indu="deprod1404",employ="delama1317",activity="desurv0048"),
              US_YC = list(M03 = "us03mgov", Y02 = "us02ygov", Y03 = "us03ygov", Y05 = "us05ygov", Y07 = "us07ygov", Y10 = "us10ygov", Y30 = "us30ygov", irs10Y = "usd10yswap", irs7Y = "usd7yswap",
                           irs5Y = "usd7yswap", irs3Y = "usd2yswap", irs2Y = "usd2yswap", Libor3M = "us3mlibor", Y10_infl_exp = "us10ybei", Y5_infl_exp = "usbkeveny05",
                           Y5Y5_infl_exp = "usbkeven5f5", MortagageRate_30Y = "usrate0218"),
              DE_YC = list(M03 = "de3mgov", Y01 = "de1ygov", Y02 = "de2ygov", Y03 = "de3ygov", Y05 = "de5ygov", Y07 = "de7ygov", Y10 = "de10ygov"),
              EA_YC = list(Euribor3M = "de3mibr", irs10Y = "ih:mb:com:eur10yswap_ext"),
              CZ_YC = list(Repo_rate = "czrate0001", Tbill_3M = "imfm935fitb_pa", Pribor3M = "cz3mpribor",irs10y="ih:mb:com:cz10yswap", Y2 = "cz2ygov", Y5 = "cz5ygov", Y10 = "cz10ygov"),
              HU_YC = list(Bubor3M="hu3mibr"),
              PL_YC = list(Wibor3M = "pl3mibr", irs10Y = "pln10yswap", Y10 = "pl10ygov"),
              HU = list(gdp = "hunaac0063", GDPrealchainedSA = "clv10mnacscab1gqhu1g", PMI = "husurv0378", PMI_New_Orders = "husurv0382", retailSales_ExAutoSA = "tovvg47scai15hutu", retaiSales_AutoSA = "tovvg47scai10hutu", IndustryProdSA = "huprod0087", ConstructionSA = "hucons0089", CarRegistrationSA_ECB = "ecb_stsmhuycregpc00003abs", unRateSA = "ecb_00143241", hicpcore="i15totxnrgfoodhudx", hicp="i05cp00hudx", unemp="ecb_00143241", eurhuf="hufecbfix", EU_funds_transfers = "mioeurkas1s1creeuihu6q"),
              CH = list(ImportUSD = "cntrad1002", import_asia ="asiatrad0001",us_chinaexp="ustrad1256",us_chinaexp_mach="ustrad_7_5700_exp",eu_chinaexp="eutrad0867",china_freight="cntran0002",china_gdp="cnnaac1909"),
              CZ_labor = list(total="czinea0019",real_estate="czinea0032",prof_service="czinea0032",manufactur="czinea0023",mining="czinea0022",
                              constr="czinea0026",education="czinea0036",financial="czinea0031",health="czinea0037",information="czinea0030",trade="czinea0027",
                              transport="czinea0028", water_supply="czinea0025",electricity="czinea0024",accomodation="czinea0029",administrative="czinea0034",
                              agriculture="czinea0020",arts="czinea0038"),
              BE = list(gdp = "benaac0274", hicp="bepric1354", unem="belama0368"),
              TR = list(GDPrealSA = "trnaac0556", GDPrealSA_USD = "trnygdpmktpsakdgemquar", GDPcurrent = "trnaac0405", GDPcurrentSA = "trnygdpmktpsacngemquar", GDPcurrentSA_USD = "trnygdpmktpsacdgemquar", ImportsRealSA = "trnaac0555", ImportsRealUSD_SA = "trnaac0555", CPI = "trpric0057", CA_balance = "trbopa0001", CA_InvestIncomeBalUSD = "trbopa0044", REER_EU = "trecexri0008", REER_CBTR = "trexri0001", REER_BIS = "bis_mrbtr", REER_JPM = "tu9rexbc10", InterbankRatesON = "trrate0043"),
              PRICE=list(cp_usa="imfq111pcpi_ix",cp_can="imfq156pcpi_ix",cp_aus="imfq193pcpi_ix",cp_mex="imfq273pcpi_ix",
                         cp_de="imfq134pcpi_ix",cp_fr="imfq132pcpi_ix",
                         cp_it="imfq136pcpi_ix", cp_es="imfq184pcpi_ix",cp_uk="imfq112pcpi_ix",cp_jp="imfq158pcpi_ix",
                         cp_kor="imfq542pcpi_ix",cp_chil="imfq228pcpi_ix",cp_indon="imfq536pcpi_ix",
                         cp_ch="imfq924pcpi_ix",cp_br="imfq223pcpi_ix",cp_ind="imfq534pcpi_ix",cp_ru="imfq922pcpi_ix",
                         cp_tr="imfq186pcpi_ix",cp_souta="imfq199pcpi_ix",cp_cz="imfq935pcpi_ix",cp_pl="imfq964pcpi_ix",
                         cp_hu="imfq944pcpi_ix",cp_be="imfq124pcpi_ix"),
              BIS = list(cred_advance = "bis_q5rpamusda", cred_emer="bis_q4tpamusda",cred_usa="bis_quspamusda",
                         cred_can="bis_qcapamxdca",cred_aus="bis_qaupamxdca",cred_mex="bis_qmxpamxdca",cred_arg="bis_qarpamxdca",
                         cred_de="bis_qdepamxdca",cred_fr="bis_qfrpamxdca", cred_it="bis_qitpamxdca",
                         cred_es="bis_qespamxdca", cred_uk="bis_qgbpamxdca",cred_jp="bis_qjppamxdca",
                         cred_kor="bis_qkrpamxdca",cred_chil="bis_qclpamxdca",cred_indon="bis_qidpamxdca",
                         cred_ch="bis_qcnpamxdca", cred_br="bis_qbrpamxdca",cred_ind="bis_qinpamxdca",
                         cred_ru="bis_qrupamxdca",cred_tr="bis_qtrpamxdca",cred_souta="bis_qzapamxdca",
                         cred_cz="bis_qczpamxdca",cred_pl="bis_qplpamxdca",cred_hu="bis_qhupamxdca",
                         cred_be="bis_qbepamxdca",res_usa="bis_qusn628",
                         res_can="bis_qcan628",res_aus="bis_qaun628",res_mex="bis_qmxn628",
                         res_de="bis_qden628",res_fr="bis_qfrn628",
                         res_it="bis_qitn628",res_es="bis_qesn628",res_uk="bis_qgbn628",res_jp="bis_qjpn628",
                         res_kor="bis_qkrn628",res_chil="bis_qcln628",res_indon="bis_qidn628",
                         res_ch="bis_qcnn628",res_br="bis_qbrn628",res_ind="bis_qinn628",res_ru="bis_qrun628",
                         res_tr="bis_qtrn628",res_souta="zabisrpp0001",res_cz="bis_qczn628",res_pl="bis_qpln628",
                         res_hu="bis_qhun628",res_be="bebisrpp0001",crgd_usa="bis_quspaa",
                         crgd_can="bis_qcapaa",crgd_aus="bis_qaupaa",crgd_mex="bis_qmxpaa",
                         crgd_de="bis_qdepaa",
                         crgd_fr="bis_qfrpaa",crgd_it="bis_qitpaa",crgd_es="bis_qespaa",crgd_uk="bis_qgbpaa",
                         crgd_jp="bis_qjppaa",crgd_kor="bis_qkrpaa",crgd_chil="bis_qclpaa",crgd_indon="bis_qidpaa",
                         crgd_ch="bis_qcnpaa",crgd_br="bis_qbrpaa",crgd_ind="bis_qinpaa",
                         crgd_ru="bis_qrupaa",crgd_tr="bis_qtrpaa",crgd_souta="bis_qzapaa",crgd_cz="bis_qczpaa",
                         crgd_pl="bis_qplpaa",crgd_hu="bis_qhupaa",crgd_be="bis_qbepaa",cpi_usa="bis_mus771",
                         cpi_de="bis_mde771",cpi_fr="bis_mfr771",cpi_it="bis_mit771",cpi_es="bis_mes771",
                         cpi_uk="bis_mgb771",cpi_jp="bis_mjp771",cpi_ch="bis_mcn771",cpi_br="bis_mbr771",
                         cpi_ind="bis_min771", cpi_ru="bis_mru771",cpi_tr="bis_mtr771",cpi_souta="bis_mza771",
                         cpi_cz="bis_mcz771",cpi_pl="bis_mpl771",cpi_hu="bis_mhu771",cpi_be="bis_mbe771"),
              EQUITY = list(spx="sp500_500",spx_pe_ca="sp500_shpe",spx_pe_fw="sp500_12mfwpe",spx_pe_now="sp500_mpe",
                            spx_sales="sp500_sales"),
              WORLD=list(global_industry="worldiptotsakdgemmonth",global_trade="worldtrad0001"))
  
  
  dat_raw <- lapply(tic, function(x) FetchTimeSeries(unlist(x)))
  
  na <- rapply(dat_raw, getIsError)
  na[na]
  
  # test missing data
  # for(i in seq_along(dat_raw)) for(j in seq_along(dat_raw[[i]])) { message(names(dat_raw[[i]][j])); as.xts(dat_raw[[i]][[j]]) }
  
  dat <- lapply(dat_raw, function(x) lapply(x, as.xts))
  dat <- lapply(dat, function(x) do.call("merge.xts", x))
  
  #Data manipulation
  dat1 <- dat
  tic <- tic
  for (i in seq_along(dat1)) colnames(dat1[[i]]) <- paste0(names(tic)[i], "_", names(tic[[i]]))
  
  #Data
  dt_X <- merge.xts(fapply_aggregate(dat1$DE[,"DE_IFO_MANU"], f = minf, per = "quarter", eop = F), # GDP
                    fapply_aggregate(dat1$DE[,"DE_IFO_CONST"], f = minf, per = "quarter", eop = F), # GDP
                    fapply_aggregate(dat1$DE[,"DE_IFO_SERV"], f = minf, per = "quarter", eop = F), # GDP
                    fapply_aggregate(dat1$DE[,"DE_IFO_TRADE"], f = minf, per = "quarter", eop = F), # GDP
    fapply_aggregate(dat1$DE[,"DE_IFO_DE"], f = minf, per = "quarter", eop = F),
    100*(diff.xts(log(fapply_aggregate(dat1$DE[,"DE_vacancies_DE"], f = minf, per = "quarter", eop = F)),lag = 1)),
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_czcorploan"], f = minf, per = "quarter", eop = F))))), # corploan
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_ind_prod"], f = minf, per = "quarter", eop = F))))), # industrial production
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_czvacan"], f = minf, per = "quarter", eop = F))))), # vacations
    #100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_empl_rate"], f = minf, per = "quarter", eop = T))))), # empl
    diff.xts(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_czunem"], f = minf, per = "quarter", eop = F))), # unem
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_new_car"], f = minf, per = "quarter", eop = F))))), # new car
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_retail"], f = minf, per = "quarter", eop = F))))), # retail
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_bus_surv"], f = minf, per = "quarter", eop = F))))), # trade balance
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_empl_m"], f = minf, per = "quarter", eop = F))))), # employment
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_cons_credit"], f = minf, per = "quarter", eop = F))))), # consumer credit
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_hous_lend"], f = minf, per = "quarter", eop = F))))), # house lending
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_machinery"], f = minf, per = "quarter", eop = F))))), # machinery
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_car_prod"], f = minf, per = "quarter", eop = F))))), # car production
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_steel_prod"], f = minf, per = "quarter", eop = F))))), # steeel production
    #100*(diff.xts((fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_ec_recentprod"], f = minf, per = "quarter", eop = T))))), # EC
    #100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_euro_induconf"], f = minf, per = "quarter", eop = T))))), # EC
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_ec_ecsent"], f = minf, per = "quarter", eop = F))))), # EC
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_czso_ecsent"], f = minf, per = "quarter", eop = F))))), # EC
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_transport"], f = minf, per = "quarter", eop = F))))), # transport
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_office_sup"], f = minf, per = "quarter", eop = F))))), # office
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_car_balance"], f = minf, per = "quarter", eop = F))))), # balance
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_machinery_exp"], f = minf, per = "quarter", eop = F))))), # machinery
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_cz_warehouse"], f = minf, per = "quarter", eop = F))))), # warehouse
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_cz_freight"], f = minf, per = "quarter", eop = F))))), # freight
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_electricity_cons"], f = minf, per = "quarter", eop = F))))), # electricity
    100*(diff.xts(log(fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_cz_advertise"], f = minf, per = "quarter", eop = F))))),
    fapply_stl(fapply_aggregate(dat1$CZ[,"CZ_PMI_CZ"], f = minf, per = "quarter", eop = F)),# PMI_CZ
    fapply_stl(fapply_aggregate(dat1$EA[,"EA_PMI_Man"], f = minf, per = "quarter", eop = F)),# PMI_EMU
    fapply_stl(fapply_aggregate(dat1$EA[,"EA_PMI_Com"], f = minf, per = "quarter", eop = F))) # PMI_EMU)
  # advertise
  matur <- index(dt_X)%m+% months(2)#stejny format dat
  dt_X <- xts(as.matrix(dt_X),order.by = matur)
  
}
#merging datasets
#dataset <- merge.xts(gdp_select,retail_xts,const_xts,int_trade_xts,unem_xts,dt_X)
#dataset <- merge.xts(gdp_select,unem_xts,dt_X)
dataset <- merge.xts(gdp_select,manu_xts,retail_xts,service_xts,const_xts,unem_xts,dt_X) # pridat int_trade_xts az budou data aktual
pokus <- dataset[paste0("/",aaa),]# data do posledniho data pro odhad modelu
pokus<- rbind(pokus,xts(matrix(colMeans(dataset[paste0(bbb,"/",(as.Date(bbb) %m+% months(2))),],na.rm = T),nrow = 1),order.by = as.Date(bbb))) # zprumerovani vsech dat nad posledni datum pro odhad modelu
dataset <- pokus


saveRDS(dataset,"dataset.rds") 
#-----------------------------------------------------------------------------------#
#Initial transformation of data#-------
#-----------------------------------------------------------------------------------#
{
  
  
  dataset.ts <- xts2ts(dataset)
  dataset <- dataset["2004-06-01/",]
  
}

#-----------------------------------------------------------------------------------#
#Data for model#
#-----------------------------------------------------------------------------------#
  dataset_select  <- dataset[, !(colnames(dataset) %in% c("gdp_emu","gdp_vol"))]#selects explanatory variables
  dataset_select <- scale(dataset_select)
  colMeans(dataset_select[paste0(bbb, "/"),],na.rm = T)
  dataset_select <- dataset_select[,!is.na(colMeans(dataset_select[paste0(bbb, "/"),],na.rm = T))]# selects variables that we have at least one observation in the period starting with the date bbb (found out by mean calculation)
  #dataset_select <- dataset_select[,!is.na(dataset_select[as.Date("2009-06-01")])]
  dataset_select_0 <- na.trim(dataset_select[paste("/",aaa,sep = ""),])#selects observations withou NA and up till the last quater (for the model to be estimated from)
  dataset_select_1 <- colMeans(dataset_select[paste0(bbb, "/"),],na.rm = T)#selects observations for the current quater (input for the forecast)
  dataset_expl_0 <- dataset[index(dataset_select_0), "gdp_emu"]#selects explanatory variable (for the model to be estimated from)
  
  
  x <- dataset_select_0 
  y <- dataset_expl_0
  x1 <- coredata(dataset_select_1) 
  x <- coredata(x)
  y <- coredata(y)
  
  # Select most important data
  #map_soub_emu <- readRDS(file = "map_soub_emu.rds")
  #x <- x[,map_soub_emu[,1]]
  #x1 <- x1[map_soub_emu[,1]]
  
  #Select soft data
  #map_soub <- matrix(c(c("EA_PMI_Com","DE_IFO_MANU","DE_IFO_CONST","DE_IFO_SERV","DE_IFO_TRADE","DE_IFO_DE","CZ_PMI_CZ","CZ_bus_surv"),c(c("EMU PMI","IFO_MANU","IFO_CONST","IFO_SERV","IFO_TRADE","IFO"),c("PMI CZ","CZSO Sentiment"))),ncol = 2)
  #x <- x[,map_soub[,1]]
  #x1 <- x1[map_soub[,1]]
  
#-----------------------------------------------------------------------------------#
## otestujeme konkurencni 'bsts' modely a vratime chyby predpovedi na jedno obdobi dopredu-----
#-----------------------------------------------------------------------------------#

iter <- 10e3
burn <- floor(iter*0.5)
{
  ## PRIKLAD 1
  # obecne definujeme zakladni 'bsts' model
  
  
  
  BSTS_models <- BSTS_TEST_CONSTRUCTOR(y = y, x = x,
                                       state_specification = AddLocalLevel(state.specification = list(), y = y),
                                       n_iterations = iter, n_burn = burn)
  BSTS_models_1 <- BSTS_TEST_CONSTRUCTOR(y = y, x = x,
                                         state_specification = AddAr(state.specification = list(), y = y),
                                         n_iterations = iter, n_burn = burn)
  
  # odhadneme jednotlive konkretni modely a ulozime chyby predpovedi
  BSTS_errors <- lapply(1:7, BSTS_models)
  BSTS_errors <- lapply(BSTS_errors, "[[", "errors")
  BSTS_errors <- lapply(BSTS_errors, "[", ,-1)
  
  BSTS_errors1 <- lapply(1:7, BSTS_models_1)
  BSTS_errors1 <- lapply(BSTS_errors1, "[[", "errors")
  BSTS_errors1 <- lapply(BSTS_errors1, "[", ,-1)
  
  BSTS_errors <- c(BSTS_errors, BSTS_errors1)
  
  # spocteme statistiky chyb predpovedi
  BSTS_result <- lapply(BSTS_errors, BSTS_MODSELECT_CRITERIA)
  
  # porovnani modelu
  BSTS_models_comparison <- BSTS_MODSELECT_SUMMARY(BSTS_result, show_model = NA)
  
  # zvyraznime model cislo 1
  BSTS_models_comparison <- BSTS_MODSELECT_SUMMARY(BSTS_result, show_model = 10)
}
#-----------------------------------------------------------------------------------#
## Odhadneme nejlepsi model------------
#-----------------------------------------------------------------------------------#
#Adding constant
#x <- cbind(x,rep(1,nrow(x)))
#colnames(x)[ncol(x)] <- "constant"
#x1 <- c(x1,1)
#names(x1)[length(x1)] <- "constant"
iter <- 10e3
burn <- floor(iter*0.5)
ss <- list()
ss <- AddLocalLevel(state.specification = ss, y = y)
#ss <- AddAr(state.specification = list(), y = y)
pr <- SpikeSlabPrior(x, y, expected.model.size = 6)



# odhadneme nejlepsi model, v tomto pripade je to model s 5 vysvetlujicimi promennymi
mod.final <- bsts(formula = y ~ x-1, 
                  state.specification = ss,
                  prior = pr,
                  niter = iter)
median(apply(mod.final$one.step.prediction.errors,MARGIN = 2,median)) #median error
median(apply(abs(mod.final$one.step.prediction.errors),MARGIN = 2,median)) # median absolute error



# grafy a ostatni vystupy
plot(mod.final, "components", burn = burn)
plot(mod.final, "size", burn = burn)
plot(mod.final, "forecast.distribution", burn = burn)
plot(mod.final, "prediction.errors", burn = burn)

plot(mod.final, "predictors", burn = burn, inclusion.threshold = 0.01)
plot(mod.final, "coefficients", burn = burn)


#-----------------------------------------------------------------------------------#
# Predictions +  estimate graph----------
#-----------------------------------------------------------------------------------#

output <- predict.bsts(object = mod.final,newdata =matrix(x1,nrow = 1),burn = burn) #prediction


#Save the nowcast sumarry
#nowcast_sumary <- matrix(ncol=3,nrow=)
#colnames(nowcast_sumary) <- c("EMU","DE","CZ")
#saveRDS(object = nowcast_sumary,file = "P:/Transfer/Dealing/AFT/HONZABURES/nowcast_sumary.rds")
nowcast_sumary <- readRDS("P:/Transfer/Dealing/AFT/HONZABURES/nowcast_sumary.rds")
nowcast_sumary[,3] <- as.numeric(output$median)
attr(nowcast_sumary,"date") <- Sys.Date()
saveRDS(object = nowcast_sumary,file = "P:/Transfer/Dealing/AFT/HONZABURES/nowcast_sumary.rds")

#Save the nowcast developmentn
# develop_cz <- xts(x = matrix(c(as.numeric(output$interval[1,]),as.numeric(output$median),as.numeric(output$interval[2,])),nrow = 1),order.by = now())
#saveRDS(object = develop_cz,file = "develop_cz.rds")

develop_cz <- readRDS(file = "develop_cz.rds")
develop_cz <- rbind(develop_cz,xts(x = matrix(c(as.numeric(output$interval[1,]),as.numeric(output$median),as.numeric(output$interval[2,])),nrow = 1),order.by = now()))
colnames(develop_cz) <- c("2.5%","Median","97.5%")
saveRDS(object = develop_cz,file = "develop_cz.rds")

#Plot of Development

png(filename = "u:/Dokumenty/Rdata/GDPnow_emu/czdevelopment.png", width = 600, height = 750)
op <- par
col2rgb(CSOB_colors(2))
mycol <- rgb(0, 51, 102, max = 255, alpha = 30, names = "blue50")
par(cex=1.4)
plot(x = index(develop_cz),y = develop_cz$Median, type = "b",pch=16,bty="n",col=CSOB_colors(1),ylim=c(-0.9,0.9),ylab = "qoqa %",xlab = "", main = "" )
polygon(x = c(index(develop_cz),index(develop_cz)[length(index(develop_cz)):1]),
        y = c(as.numeric(develop_cz$`2.5%`),as.numeric(develop_cz$`97.5%`)[length(index(develop_cz)):1]),col =mycol, border = NA )
text(x = index(develop_cz),y =(develop_cz$Median-0.1),labels = round(as.numeric(develop_cz$Median),digits = 2))
abline(h = 0,col="red",lty=2)
par(op)
dev.off()

#plot of comparison with past
temp_date <- "2010-01-01/"
plot(y = as.numeric(dataset[temp_date,"gdp_emu"]),x = index(dataset[temp_date,]),bty="null",ylab = "Mezikvartalni rust HDP CZ",xlab = "")
abline(h = as.numeric(output$median),col="grey40",lty=2)
text(paste("Odhad mezikvartalniho rustu 2019 Q3:",round(output$median,digits = 2),"%"),y =output$median-0.5,x = index(dataset)[62],col ="grey40",adj = 1)
points(x = last(index(dataset)),y =output$median,pch=19)
points(x = last(index(dataset)),y =output$median,pch=2)
#Probability above threshold

prob <- ecdf(output$distribution)(0) * 100 # funkce ecdf vrati funkci kumulativni distribuce, ktera nasledne vyhodi vysledek pro zadanou hodnotu ....0 znamena mensi pravdepodobnost, ze hodnota bude mensi nebo rovna 0  

op <- par
par(cex=1.4)
hist(output$distribution,xlab = "", ylab="Hustota pravdepodobnosti",freq = 0,main = "Nowcast odhad HDP v CZ 4Q 2022",col = "grey90",border = "white")
text(paste("Odhad mezikvartalniho rustu(median):",round(output$median,digits = 1),"%"),x =output$median+2,y = 0.1,col ="grey40",adj = 1,cex = 0.7)
text(paste("Pravdepodobnost poklesu HDP:",round((prob),digits = 1),"%"),x =output$median+2,y = 0.06,col ="grey40",adj = 1,cex = 0.7 )
#text(paste("Odhad CNB:","-9,4%"),x =output$median-2,y = 0.05,col ="grey40",adj = 1,cex = 0.7 )
abline(v=0,col="grey25",lty=2 )
mtext(text = "Zdroj: CSOB/Patria", side = 1, line = 2, adj = 1, font = 2)
par(op)

#twitter size
op <- par
par(cex=1.4)
hist(output$distribution,xlab = "", ylab="Hustota pravdepodobnosti",freq = 0,main = "Nowcast odhad HDP v CZ 2Q 2020",col = "grey90",border = "white",cex.main=0.7,cex.axis=0.7)
text(paste("Odhad rustu(median,qoq):",round(output$median,digits = 1),"%"),x =output$median-3,y = 0.07,col ="grey40",adj = 1,cex = 0.6)
text(paste("Pravdepodobnost poklesu HDP:",round((prob),digits = 1),"%"),x =output$median-3,y = 0.05,col ="grey40",adj = 1,cex = 0.6 )
text(paste("Odhad CNB:","-9,4%"),x =output$median-3,y = 0.03,col ="grey40",adj = 1,cex = 0.6 )
abline(v=0,col="grey25",lty=2 )
mtext(text = "Zdroj: CSOB/Patria", side = 1, line = 2, adj = 1, font = 2,cex = 0.5)
par(op)

png(filename = "u:/Dokumenty/Rdata/GDPnow_emu/cznowcast.png", width = 600, height = 750)
op <- par
par(cex=1.4)
hist(output$distribution,xlab = "", ylab="Probability",freq = 0,main = "Nowcast Estimate GDP,CZ 2022Q1",col = "grey90",border = "white")
text(paste("GDP growth estimate(q/q,median):",round(output$median,digits = 2),"%"),x =output$median+5,y = 0.10,col ="grey40",adj = 1)
text(paste("Probability of GDP decline:",round((prob),digits = 1),"%"),x =output$median+5,y = 0.05,col ="grey40",adj = 1 )
abline(v=0,col="grey25",lty=2 )
mtext(text = "Source: CSOB/KBC", side = 1, line = 2, adj = 1, font = 2)
par(op)
dev.off()

#-----------------------------------------------------------------------------------#
# Compare inclusion of important variables--------------
#-----------------------------------------------------------------------------------#
#getting important vaariables about threshold 1percent and their coefficients
pok <- apply(mod.final$coefficients[-1:-burn,],MARGIN = 2,FUN = function(x){
  del <- length(x)
  poc <- length(subset(x,x != 0))
  prav <- poc/del
  return(prav)
})

pok1 <- apply(mod.final$coefficients[-1:-burn,],MARGIN = 2,FUN = function(x){
  poc <- median(subset(x,x != 0))
  return(poc)
})


pok <- sort(pok[pok > l], decreasing = T)
pok <- rbind(pok,pok1[names(pok)])

pok <- t(as.data.frame(pok))
for(i in 1:length(rownames(pok))){
  rownames(pok)[i] <- substr(x = rownames(pok)[i],start = 2,stop = 1000)
}
colnames(pok) <- c("Probability","Coefficient")
pok



map_soub_emu[!map_soub_emu[,1]%in%rownames(pok),] # Is present in map file, but not significant in current model outcome
pok[(map_soub_emu[!map_soub_emu[,1]%in%rownames(pok),][,1]),]
pokus[!rownames(pok)%in%map_soub_emu[,1],] # Is missing in the map file, although significant in current model oucome



map_soub_cz <- matrix(c(c("CZ_czvacan","C","B.D","B_C","EA_PMI_Com","exp","MIG_CAG","imp","CZ_car_balance","C23","CZ_ind_prod","CZ_PMI_CZ","MIG_ING","G45"),
                     c("Vacancies","Manufacturing","Mining and electricity.","Mining and manufacturing","EMU PMI comp.",
                     "Exports","Capital goods","Imports","Car sector balance","Concrete & cement","Industry total","CZ PMI","Intermediate goods","Retail: cars")),ncol = 2)

#Adjust the map_file  so non-significant values are deleted
#map_soub_emu <- map_soub_emu[map_soub_emu[,1]%in%pokus$V2,]


#saveRDS(object = map_soub_cz,file = "map_soub_cz.rds")
map_soub_cz <- readRDS(file = "map_soub_cz.rds")




#-----------------------------------------------------------------------------------#
# Graph of explanator y variables------------
#-----------------------------------------------------------------------------------#
#av <- apply(X = dataset_select_0[,map_soub[,1]], FUN=mean,MARGIN = 2)
#st <- apply(X = dataset_select_0[,map_soub[,1]], FUN=sd,MARGIN = 2)
hodnoc <- dataset_select_1[map_soub_cz[,1]]

#hodnoc <- (akt-av)/st
names(hodnoc) <- mapvalues(x = names(hodnoc),from =map_soub_cz[,1],to = map_soub_cz[,2] )

temp_select <- match(c("Vacancies","Manufacturing","Mining and electricity.","Exports","EMU PMI comp."),names(sort(hodnoc)))
bord_fill <- rep(NA,length(hodnoc))
bord_fill[temp_select] <- "black"

text_fill <- rep(NA,length(hodnoc))
text_fill[temp_select] <- round(sort(hodnoc)[temp_select],digits = 1 )

#names(hodnoc) <- vyb_rename
png(filename = "u:/Dokumenty/Rdata/GDPnow_emu/czexplenatory.png", width = 1500, height = 950)
k <- par(no.readonly=T)
par(mar=c(8,4,4,3),cex.axis=0.7,cex=1)
p <- barplot(sort(hodnoc),space = 1,border = bord_fill, col="lightblue1", main="Main explanatory variables, deviation from average",las=2, ylab="number of st.dev.")
text(x = p,y =-0.5,labels = text_fill,col ="grey20",font = 2 )
par(k)
dev.off()


#analgraf <- dataset.ts[,vyb]
#colnames(analgraf) <- vyb_rename
#HP_decomposition(analgraf,log=F,plot=T)


#-----------------------------------------------------------------------------------#
# copmutation of y/y growth rates, putting together y/y estimates-------------
#-----------------------------------------------------------------------------------#
macrobond_gdp <- FetchTimeSeries("clv10mnacscab1gqcz1g")
macrobond_gdp <- as.xts(macrobond_gdp$clv10mnacscab1gqcz1g)

qq_gdp <- exp(diff.xts(log(macrobond_gdp), 1)) * 100 - 100
yy_gdp <- exp(diff.xts(log(macrobond_gdp), 4)) * 100 - 100
annual_gdp <- exp(diff.xts(log(apply.yearly(macrobond_gdp,mean)),1))*100-100


grates <- c(0.1,
            0.3,
            0.9,
            0.7,
            1.2,
            1,
            1,
            0.9)#urcim tempa rustu na nasledujici roky
grates <- c(0,grates)
grates <- (grates/100)+1
grates <- cumprod(grates)

newobs_names <- seq(from=last(index(na.trim(macrobond_gdp))),by = "3 months",length.out = length(grates))
newobs <- xts(x = grates,order.by = newobs_names)
colnames(newobs) <- "rust"
newobs$gdp_emuvol <- grates*as.numeric(last(na.trim(macrobond_gdp)))
newobs <- newobs[-1,"gdp_emuvol"]
newobs <- rbind(macrobond_gdp,newobs)

qq_gdp_f <- exp(diff.xts(log(newobs), 1)) * 100 - 100
yy_gdp_f <- exp(diff.xts(log(newobs), 4)) * 100 - 100
annual_gdp_f <- exp(diff.xts(log(apply.yearly(newobs,mean)),1))*100-100

newobs_yy <- (diff(newobs,lag = 4)/lag.xts(newobs,k = 4))*100
y_newobs_yy <- (diff(y_newobs,lag = 1)/lag.xts(y_newobs,k = 1))*100

tail(newobs_yy,n = length(grates))
tail(y_newobs_yy)

#-----------------------------------------------------------------------------------#
# final graphs---------------
#-----------------------------------------------------------------------------------#
temp_dat <- "1999-01-01/"
png(filename = "u:/Dokumenty/Rdata/GDPnow_emu/czexplenatory1.png", width = 1500, height = 950)
op <- par()

par(mfrow=c(2,2),cex=1)
#Vacancies
dat_a <- dataset[temp_dat,c("CZ_czvacan","gdp_emu")]
colnames(dat_a) <- c("Vacancies","GDP")
plot(x = index(dat_a),y=dat_a[,"Vacancies"],type="b", col="grey60",bty="n",xlab="",ylab="q/q(%)",lty=2)
lines(x = index(dat_a),y=dat_a[,"GDP"], col="grey15")
points(x = last(index(dat_a)),y = last(dat_a[,"Vacancies"]),pch=10,col=CSOB_colors(2))
legend(legend = c("Vacancies","GDP"),col = c("grey60","grey15"),lty = c(2,1),x = "bottomleft",bty="n")
grid(lwd = 1,col="grey80")

#Manufacturing
dat_a <- dataset[temp_dat,c("C","gdp_emu")]
colnames(dat_a) <- c("Manufacturing","GDP")
plot(x = index(dat_a),y=dat_a[,"Manufacturing"],type="b", col="grey60",bty="n",xlab="",ylab="q/q(%)",lty=2)
lines(x = index(dat_a),y=dat_a[,"GDP"], col="grey15")
points(x = last(index(dat_a)),y = last(dat_a[,"Manufacturing"]),pch=10,col=CSOB_colors(2))
legend(legend = c("Manufacturing","GDP"),col = c("grey60","grey15"),lty = c(2,1),x = "bottomleft",bty="n")
grid(lwd = 1,col="grey80")

#Energy and mining
dat_a <- dataset[temp_dat,c("B.D","gdp_emu")]
colnames(dat_a) <- c("Energy and mining","GDP")
plot(x = index(dat_a),y=dat_a[,"Energy and mining"],type="b", col="grey60",bty="n",xlab="",ylab="q/q(%)",lty=2)
lines(x = index(dat_a),y=dat_a[,"GDP"], col="grey15")
points(x = last(index(dat_a)),y = last(dat_a[,"Energy and mining"]),pch=10,col=CSOB_colors(2))
legend(legend = c("Energy and mining","GDP"),col = c("grey60","grey15"),lty = c(2,1),x = "bottomleft",bty="n")
grid(lwd = 1,col="grey80")

#Machinery
dat_a <- dataset[temp_dat,c("F","gdp_emu")]
colnames(dat_a) <- c("Construction","GDP")
plot(x = index(dat_a),y=dat_a[,"Construction"],type="b", col="grey60",bty="n",xlab="",ylab="q/q(%)",lty=2)
lines(x = index(dat_a),y=dat_a[,"GDP"], col="grey15")
points(x = last(index(dat_a)),y = last(dat_a[,"Construction"]),pch=10,col=CSOB_colors(2))
legend(legend = c("Construction","GDP"),col = c("grey60","grey15"),lty = c(2,1),x = "bottomleft",bty="n")
grid(lwd = 1,col="grey80")
par <- op

#Machinery
dat_a <- dataset["2015-01-01/",c("H.N_STS","gdp_emu")]
colnames(dat_a) <- c("Transport and administartive services","GDP")
plot(x = index(dat_a),y=dat_a[,"Transport and administartive services"],type="b", col="grey60",bty="n",xlab="",ylab="q/q(%)",lty=2)
lines(x = index(dat_a),y=dat_a[,"GDP"], col="grey15")
points(x = last(index(dat_a)),y = last(dat_a[,"Transport and administartive services"]),pch=10,col=CSOB_colors(2))
legend(legend = c("Transport and administartive services","GDP"),col = c("grey60","grey15"),lty = c(2,1),x = "bottomleft",bty="n")
grid(lwd = 1,col="grey80")
dev.off()
