##
## Script for reading in ACE data and testing stuff.
##

# StaIClassi routines are expected to be in the parent folder
addpath('../')


# Read in ACE SWICS/SWIMS data
if (0)

  swics_files = glob("ACE_SWICS-SWIMS/ACE_SWICS_Data*.txt");

  stime_swics = [];
  nHe2 = vHe2 = vC5 = vO6 = vFe10 = vthHe2 = vthC5 = vthO6 = vthFe10 = He_q = C5_q = O6_q = Fe10_q =  [];
  C6to5 = O7to6 = avqC = avqO = avqFe = C6to5_q = O7to6_q = avqC_q = avqO_q = avqFe_q = [];
  SW_type = FetoO = FetoO_q = [];

  for i=1:numel(swics_files)

    # headers are not consistent in length, so count lines to 'DATE', close file; then textread()
    fid = fopen(swics_files{i});
    hdr_lines=0;
    do hdr_lines++; until strncmp(fgets(fid), 'BEGIN DATA', 10)
    fclose(fid);


    # convert date strings into Octave serial times
    [year_tmp, doy_tmp, hour_tmp, min_tmp, sec_tmp, fyear_tmp, fdoy_tmp, ACEepoch_tmp, ...
     nHe2_tmp, vHe2_tmp, vC5_tmp, vO6_tmp, vFe10_tmp, vthHe2_tmp, vthC5_tmp, vthO6_tmp, vthFe10_tmp, He_q_tmp, C5_q_tmp, O6_q_tmp, Fe10_q_tmp,...
     C6to5_tmp, O7to6_tmp, avqC_tmp, avqO_tmp, avqFe_tmp, C6to5_q_tmp, O7to6_q_tmp, avqC_q_tmp, avqO_q_tmp, avqFe_q_tmp, SW_type_tmp, FetoO_tmp, FetoO_q_tmp] = ...
       textread(swics_files{i}, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HEADERLINES', hdr_lines);

    # I'm not sure this is any faster than an explicit loop for each line
    stime_tmp = datenum(year_tmp, 0, doy_tmp, hour_tmp, min_tmp, round(sec_tmp));

    # append the vectors, then (re)sort by serial time
    [stime_swics, idx] = sort([stime_swics; stime_tmp]);
    nHe2 = [nHe2; nHe2_tmp](idx);
    vHe2 = [vHe2; vHe2_tmp](idx);
    vC5 = [vC5; vC5_tmp](idx);
    vO6 = [vO6; vO6_tmp](idx);
    vFe10 = [vFe10; vFe10_tmp](idx);
    vthHe2 = [vthHe2; vthHe2_tmp](idx);
    vthC5 = [vthC5; vthC5_tmp](idx);
    vthO6 = [vthO6; vthO6_tmp](idx);
    vthFe10 = [vthFe10; vthFe10_tmp](idx);
    He_q = [He_q; He_q_tmp](idx);
    C5_q = [C5_q; C5_q_tmp](idx);
    O6_q = [O6_q; O6_q_tmp](idx);
    Fe10_q = [Fe10_q; Fe10_q_tmp](idx);
    C6to5 = [C6to5; C6to5_tmp](idx);
    O7to6 = [O7to6; O7to6_tmp](idx);
    avqC = [avqC; avqC_tmp](idx);
    avqO = [avqO; avqO_tmp](idx);
    avqFe = [avqFe; avqFe_tmp](idx);
    C6to5_q = [C6to5_q; C6to5_q_tmp](idx);
    O7to6_q = [O7to6_q; O7to6_q_tmp](idx);
    avqC_q = [avqC_q; avqC_q_tmp](idx);
    avqO_q = [avqO_q; avqO_q_tmp](idx);
    avqFe_q = [avqFe_q; avqFe_q_tmp](idx);
    SW_type = [SW_type; SW_type_tmp](idx);
    FetoO = [FetoO; FetoO_tmp];
    FetoO_q = [FetoO_q; FetoO_q_tmp];



    printf("\rDone reading %s (file %d of %d)", swics_files{i}, i, numel(swics_files));
    fflush(stdout);

  endfor

  clear *_tmp i j hdr_lines idx fid


endif


# Read in ACE Magnetometer data
if (0)

  mag_files = glob("ACE_MAG/ACE_MAG_Data*.txt");

  stime_mag = [];
  Br = Bt = Bn = Bmag = Delta = Lambda = [];
  Bgse_x = Bgse_y = Bgse_z = Bgsm_x = Bgsm_y = Bgsm_z = [];
  dBrms = sigma_B = Quality = [];

  for i=1:numel(mag_files)

    # headers are not consistent in length, so count lines to 'DATE', close file; then textread()
    fid = fopen(mag_files{i});
    hdr_lines=0;
    do hdr_lines++; until strncmp(fgets(fid), 'BEGIN DATA', 10)
    fclose(fid);

    # convert date strings into Octave serial times
    [year_tmp, doy_tmp, hour_tmp, min_tmp, sec_tmp, fyear_tmp, fdoy_tmp, ACEepoch_tmp, ...
     Br_tmp, Bt_tmp, Bn_tmp, Bmag_tmp, Delta_tmp, Lambda_tmp, Bgse_x_tmp, Bgse_y_tmp, Bgse_z_tmp, Bgsm_x_tmp, Bgsm_y_tmp, Bgsm_z_tmp, dBrms_tmp, sigma_B_tmp, Quality_tmp] = ...
       textread(mag_files{i}, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HEADERLINES', hdr_lines);

    # I'm not sure this is any faster than an explicit loop for each line
    stime_tmp = datenum(year_tmp, 0, doy_tmp, hour_tmp, min_tmp, round(sec_tmp));


    # append the vectors, then (re)sort by serial time
    [stime_mag, idx] = sort([stime_mag; stime_tmp]);
    Br = [Br; Br_tmp](idx);
    Bt = [Bt; Bt_tmp](idx);
    Bn = [Bn; Bn_tmp](idx);
    Bmag = [Bmag; Bmag_tmp](idx);
    Delta = [Delta; Delta_tmp](idx);
    Lambda = [Lambda; Lambda_tmp](idx);
    Bgse_x = [Bgse_x; Bgse_x_tmp](idx);
    Bgse_y = [Bgse_y; Bgse_y_tmp](idx);
    Bgse_z = [Bgse_z; Bgse_z_tmp](idx);
    Bgsm_x = [Bgsm_x; Bgsm_x_tmp](idx);
    Bgsm_y = [Bgsm_y; Bgsm_y_tmp](idx);
    Bgsm_z = [Bgsm_z; Bgsm_z_tmp](idx);
    dBrms = [dBrms; dBrms_tmp];
    sigma_B = [sigma_B; sigma_B_tmp];
    Quality = [Quality; Quality_tmp];


    printf("\rDone reading %s (file %d of %d)", mag_files{i}, i, numel(mag_files));
    fflush(stdout);

  endfor

  clear *_tmp i j hdr_lines idx fid


endif



# Read in preliminary ACE SWEPAM data
if (0)

  swe_pre_files = glob("ACE_SWEPAM/AC_K1_SWE_*.txt");

  stime_swe_pre = [];
  speed_swe_pre = temp_swe_pre = rho_swe_pre = HeHratio_swe_pre = [];

  for i=1:numel(swe_pre_files)

    # headers are not consistent in length, so count lines to 'DATE', close file; then textread()
    fid = fopen(swe_pre_files{i});
    hdr_lines=0;
    do hdr_lines++; until strncmp(fgets(fid), 'dd-mm-yyyy', 10)
    fclose(fid);

    # convert date strings into Octave serial times
    [date_str_tmp, time_str_tmp, rho_sw_tmp, speed_sw_tmp, HeHratio_sw_tmp, temp_sw_tmp] = ...
       textread(swe_pre_files{i}, '%s %s %f %f %f %f', 'HEADERLINES', hdr_lines);

    nrecs = 6912; % this is the number of hours between 8/1/2010, 0h and 5/15/2011, 23h

    # I'm not sure this is any faster than an explicit loop for each line
    stime_tmp = datenum(datevec([cell2mat(date_str_tmp(1:6912)), ...
                                repmat(' ', 6912, 1), ...
                                cell2mat(time_str_tmp(1:6912))], "dd-mm-yyyy HH:MM:SS.FFF"));


    # append the vectors, then (re)sort by serial time
    [stime_swe_pre, idx] = sort([stime_swe_pre; stime_tmp]);
    rho_swe_pre = [rho_swe_pre; rho_sw_tmp](idx);
    speed_swe_pre = [speed_swe_pre; speed_sw_tmp](idx);
    HeHratio_swe_pre = [HeHratio_swe_pre; HeHratio_sw_tmp](idx);
    temp_swe_pre = [temp_swe_pre; temp_sw_tmp](idx);


    printf("\rDone reading %s (file %d of %d)", swe_pre_files{i}, i, numel(swe_pre_files));
    fflush(stdout);

  endfor

  clear *_tmp i j hdr_lines idx fid

endif



# Read in level 2 ACE SWEPAM data
if (0)

  swe_files = glob("ACE_SWEPAM/ACE_swepam_level2*.txt");

  stime_swe = [];
  speed_swe = temp_swe = rho_swe = HeHratio_swe = [];

  for i=1:numel(swe_files)

    # headers are not consistent in length, so count lines to 'DATE', close file; then textread()
    fid = fopen(swe_files{i});
    hdr_lines=0;
    do hdr_lines++; until strncmp(fgets(fid), 'BEGIN DATA', 10)
    fclose(fid);

    # convert date strings into Octave serial times
    [year_tmp, doy_tmp, hour_tmp, min_tmp, sec_tmp, fyear_tmp, fdoy_tmp, ACEepoch_tmp, ...
     rho_sw_tmp, temp_sw_tmp, HeHratio_sw_tmp, speed_sw_tmp] = ...
       textread(swe_files{i}, '%f %f %f %f %f %f %f %f %f %f %f %f', 'HEADERLINES', hdr_lines);

    # I'm not sure this is any faster than an explicit loop for each line
    stime_tmp = datenum(year_tmp, 0, doy_tmp, hour_tmp, min_tmp, round(sec_tmp));


    # append the vectors, then (re)sort by serial time
    [stime_swe, idx] = sort([stime_swe; stime_tmp]);
    rho_swe = [rho_swe; rho_sw_tmp](idx);
    temp_swe = [temp_swe; temp_sw_tmp](idx);
    HeHratio_swe = [HeHratio_swe; HeHratio_sw_tmp](idx);
    speed_swe = [speed_swe; speed_sw_tmp](idx);


    printf("\rDone reading %s (file %d of %d)", swe_files{i}, i, numel(swe_files));
    fflush(stdout);

  endfor

  clear *_tmp i j hdr_lines idx fid


endif


# Merge preliminary and L2 SWEPAM data to obtain best coverage
if (0)

  pre_start = datenum(2010,8,1);

  stime_swe_merged = [stime_swe(stime_swe<pre_start); stime_swe_pre(stime_swe_pre>=pre_start)];
  rho_swe_merged = [rho_swe(stime_swe<pre_start); rho_swe_pre(stime_swe_pre>=pre_start)];
  temp_swe_merged = [temp_swe(stime_swe<pre_start); temp_swe_pre(stime_swe_pre>=pre_start)];
  HeHratio_swe_merged = [HeHratio_swe(stime_swe<pre_start); HeHratio_swe_pre(stime_swe_pre>=pre_start)];
  speed_swe_merged = [speed_swe(stime_swe<pre_start); speed_swe_pre(stime_swe_pre>=pre_start)];

  clear pre_start


endif


#
# Read in Alysha Reinard's ACE L2 data from ~1998 to ~2009
#

# Read in SWEPAM data
if (0)

  swe_files = glob("ACE_SWEPAM/ACE_swepam_l2_1998to2009.txt");

  stime_swe = [];
  Hspeed = Htemp = Hrho = HeHratio = [];
  vx_gse = vy_gse = vz_gse = vr = vt = vn = vx_gsm = vy_gsm = vz_gsm = [];
  rx_gse_swe = ry_gse_swe = rz_gse_swe = rx_gsm_swe = ry_gsm_swe = rz_gsm_swe = [];

  for i=1:numel(swe_files)

    # headers are not consistent in length, so count lines to 'DATE', close file; then textread()
    fid = fopen(swe_files{i});
    hdr_lines=0
    #do hdr_lines++; until strncmp(fgets(fid), 'BEGIN DATA', 10)
    fclose(fid);
    fflush(stdout);


    # convert date strings into Octave serial times
    [year_tmp, doy_tmp, hour_tmp, min_tmp, sec_tmp, fyear_tmp, fdoy_tmp, ACEepoch_tmp, ...
     Hrho_tmp, Htemp_tmp, HeHratio_tmp, Hspeed_tmp, ...
     vx_gse_tmp, vy_gse_tmp, vz_gse_tmp, vr_tmp, vt_tmp, vn_tmp, vx_gsm_tmp, vy_gsm_tmp, vz_gsm_tmp, ...
     rx_gse_tmp, ry_gse_tmp, rz_gse_tmp, rx_gsm_tmp, ry_gsm_tmp, rz_gsm_tmp, ...
     c28_tmp, c29_tmp, c30_tmp, c31_tmp, c32_tmp, c33_tmp, c34_tmp, c35_tmp, c36_tmp] = ...
       textread(swe_files{i}, ["%f %f %f %f %f %f %f %f %f %f %f %f ",...
                               "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
                               "%f %f %f %f %f %f %f %f %f"], 'HEADERLINES', hdr_lines);

    # I'm not sure this is any faster than an explicit loop for each line
    stime_tmp = datenum(year_tmp, 0, doy_tmp, hour_tmp, min_tmp, round(sec_tmp));


    # append the vectors, then (re)sort by serial time
    [stime_swe, idx] = sort([stime_swe; stime_tmp]);
    Hrho = [Hrho; Hrho_tmp](idx);
    Htemp = [Htemp; Htemp_tmp](idx);
    HeHratio = [HeHratio; HeHratio_tmp](idx);
    Hspeed = [Hspeed; Hspeed_tmp](idx);
    vx_gse = [vx_gse; vx_gse_tmp](idx);
    vy_gse = [vy_gse; vy_gse_tmp](idx);
    vz_gse = [vz_gse; vz_gse_tmp](idx);
    vr = [vr; vr_tmp](idx);
    vt = [vt; vt_tmp](idx);
    vn = [vn; vn_tmp](idx);
    vx_gsm = [vx_gsm; vx_gsm_tmp](idx);
    vy_gsm = [vy_gsm; vy_gsm_tmp](idx);
    vz_gsm = [vz_gsm; vz_gsm_tmp](idx);
    rx_gse_swe = [rx_gse_swe; rx_gse_tmp](idx);
    ry_gse_swe = [ry_gse_swe; ry_gse_tmp](idx);
    rz_gse_swe = [rz_gse_swe; rz_gse_tmp](idx);
    rx_gsm_swe = [rx_gsm_swe; rx_gsm_tmp](idx);
    ry_gsm_swe = [ry_gsm_swe; ry_gsm_tmp](idx);
    rz_gsm_swe = [rz_gsm_swe; rz_gsm_tmp](idx);


    printf("\rDone reading %s (file %d of %d)", swe_files{i}, i, numel(swe_files));
    fflush(stdout);

  endfor

  # convert bad measurements to NaNs for subsequent processing
  Hrho(Hrho == -9999.9) = nan;
  HeHratio(HeHratio == -9999.9) = nan;
  Htemp(Htemp == -9999.9) = nan;
  Hspeed(Hspeed == -9999.9) = nan;
  vx_gse(vx_gse == -9999.9) = nan;
  vy_gse(vy_gse == -9999.9) = nan;
  vz_gse(vz_gse == -9999.9) = nan;
  vr(vr == -9999.9) = nan;
  vt(vt == -9999.9) = nan;
  vn(vn == -9999.9) = nan;
  vx_gsm(vx_gsm == -9999.9) = nan;
  vy_gsm(vy_gsm == -9999.9) = nan;
  vz_gsm(vz_gsm == -9999.9) = nan;


  printf("\n");
  clear *_tmp i j hdr_lines idx fid


endif



# Read in MAG data
if (0)

  mag_files = glob("ACE_MAG/ACE_mag_l2_1997to2009.txt");

  stime_mag = [];
  Br = Bt = Bn = Bmag = Delta = Lambda = [];
  Bx_gse = By_gse = Bz_gse = Bx_gsm = By_gsm = Bz_gsm = [];
  dBrms = sigma_B = frac_good = nvec = Quality = [];
  rx_gse_mag = ry_gse_mag = rz_gse_mag = rx_gsm_mag = ry_gsm_mag = rz_gsm_mag = [];

  for i=1:numel(mag_files)

    # headers are not consistent in length, so count lines to 'DATE', close file; then textread()
    fid = fopen(mag_files{i});
    hdr_lines=0;
    do hdr_lines++; until strncmp(fgets(fid), 'BEGIN DATA', 10)
    fclose(fid);


    # convert date strings into Octave serial times
    [year_tmp, doy_tmp, hour_tmp, min_tmp, sec_tmp, fyear_tmp, fdoy_tmp, ACEepoch_tmp, SCclock, ...
     Br_tmp, Bt_tmp, Bn_tmp, Bmag_tmp, Delta_tmp, Lambda_tmp, ...
     Bx_gse_tmp, By_gse_tmp, Bz_gse_tmp, Bx_gsm_tmp, By_gsm_tmp, Bz_gsm_tmp, ...
     dBrms_tmp, sigma_B_tmp, frac_good_tmp, nvec_tmp, Quality_tmp, ...
     rx_gse_tmp, ry_gse_tmp, rz_gse_tmp, rx_gsm_tmp, ry_gsm_tmp, rz_gsm_tmp] = ...
       textread(mag_files{i}, ["%f %f %f %f %f %f %f %f %f %f ",...
                               "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
                               "%f %f %f %f %f"], 'HEADERLINES', hdr_lines);

    # I'm not sure this is any faster than an explicit loop for each line
    stime_tmp = datenum(year_tmp, 0, doy_tmp, hour_tmp, min_tmp, round(sec_tmp));


    # append the vectors, then (re)sort by serial time
    [stime_mag, idx] = sort([stime_mag; stime_tmp]);
    Br = [Br; Br_tmp](idx);
    Bt = [Bt; Bt_tmp](idx);
    Bn = [Bn; Bn_tmp](idx);
    Bmag = [Bmag; Bmag_tmp](idx);
    Delta = [Delta; Delta_tmp](idx);
    Lambda = [Lambda; Lambda_tmp](idx);
    Bx_gse = [Bx_gse; Bx_gse_tmp](idx);
    By_gse = [By_gse; By_gse_tmp](idx);
    Bz_gse = [Bz_gse; Bz_gse_tmp](idx);
    Bx_gsm = [Bx_gsm; Bx_gsm_tmp](idx);
    By_gsm = [By_gsm; By_gsm_tmp](idx);
    Bz_gsm = [Bz_gsm; Bz_gsm_tmp](idx);
    dBrms = [dBrms; dBrms_tmp](idx);
    sigma_B = [sigma_B; sigma_B_tmp](idx);
    frac_good = [frac_good; frac_good_tmp](idx);
    nvec = [nvec; nvec_tmp](idx);
    Quality = [Quality; Quality_tmp](idx);
    rx_gse_mag = [rx_gse_mag; rx_gse_tmp](idx);
    ry_gse_mag = [ry_gse_mag; ry_gse_tmp](idx);
    rz_gse_mag = [rz_gse_mag; rz_gse_tmp](idx);
    rx_gsm_mag = [rx_gsm_mag; rx_gsm_tmp](idx);
    ry_gsm_mag = [ry_gsm_mag; ry_gsm_tmp](idx);
    rz_gsm_mag = [rz_gsm_mag; rz_gsm_tmp](idx);


    printf("\rDone reading %s (file %d of %d)", mag_files{i}, i, numel(mag_files));
    fflush(stdout);

  endfor


  # convert bad measurements to NaNs for subsequent processing
  Br(Br == -999.9) = nan;
  Bt(Bt == -999.9) = nan;
  Bn(Bn == -999.9) = nan;
  Bmag(Bmag == -999.9) = nan;
  Delta(Delta == -999.9) = nan;
  Lambda(Lambda == -999.9) = nan;
  Bx_gse(Bx_gse == -999.9) = nan;
  By_gse(By_gse == -999.9) = nan;
  Bz_gse(Bz_gse == -999.9) = nan;
  Bx_gsm(Bx_gsm == -999.9) = nan;
  By_gsm(By_gsm == -999.9) = nan;
  Bz_gsm(Bz_gsm == -999.9) = nan;
  dBrms(dBrms == -999.9) = nan;
  sigma_B(sigma_B == -999.9) = nan;
  frac_good(frac_good == -999.9) = nan;
  nvec(nvec == -999.9) = nan;
  Quality(Quality == -999.9) = nan;




  printf("\n");
  clear *_tmp i j hdr_lines idx fid



endif





# Read in SWICS-SWIM data
if (0)

  swics_files = glob("ACE_SWICS-SWIMS/ACE_swics_l2_1998to2009.txt");

  stime_swics = [];
  nHe2 = vHe2 = vC5 = vO6 = vFe10 = vthHe2 = vthC5 = vthO6 = vthFe10 = [];
  He_qual = C5_qual = O6_qual = Fe10_qual = [];
  C6to5 = O7to6 = avqC = avqO = avqFe = [];
  C6to5_qual = O7to6_qual = avqC_qual = avqO_qual = avqFe_qual = SW_type = FetoO = FetoO_qual = [];



  for i=1:numel(swics_files)

    # headers are not consistent in length, so count lines to 'DATE', close file; then textread()
    fid = fopen(swics_files{i});
    hdr_lines=0;
    do hdr_lines++; until strncmp(fgets(fid), 'BEGIN DATA', 10)
    fclose(fid);


    # convert date strings into Octave serial times
    [year_tmp, doy_tmp, hour_tmp, min_tmp, sec_tmp, fyear_tmp, fdoy_tmp, ACEepoch_tmp, ...
     nHe2_tmp, vHe2_tmp, vC5_tmp, vO6_tmp, vFe10_tmp, ...
     vthHe2_tmp, vthC5_tmp, vthO6_tmp, vthFe10_tmp, ...
     He_qual_tmp, C5_qual_tmp, O6_qual_tmp, Fe10_qual_tmp, ...
     C6to5_tmp, O7to6_tmp, avqC_tmp, avqO_tmp, avqFe_tmp, ...
     C6to5_qual_tmp, O7to6_qual_tmp, avqC_qual_tmp, avqO_qual_tmp, avqFe_qual_tmp, ...
     SW_type_tmp, FetoO_tmp, FetoO_qual_tmp] = ...
       textread(swics_files{i}, ["%f %f %f %f %f %f %f %f ",...
                                 "%f %f %f %f %f %f %f %f %f %f %f %f %f ",...
                                 "%f %f %f %f %f %f %f %f %f %f %f %f %f"], 'HEADERLINES', hdr_lines);


    # SWICS time stamps are completely nonsensical...after rounding to nearest hour (error source #1), we
    # drop duplicate time stamps (error source #2), then fill missing steps with -9999.9 (error source #3)
    #stime_tmp = datenum(year_tmp, 0, doy_tmp, hour_tmp, min_tmp, round(sec_tmp));
    hour_tmp = round(hour_tmp + min_tmp/60 + sec_tmp/3600);
    stime_tmp = datenum(year_tmp, 0, doy_tmp, hour_tmp);

    stime_dups = [diff(stime_tmp) == 0; false]; # this removes 1st instance of repeated time stamp
    stime_tmp(stime_dups) = [];
    nHe2_tmp(stime_dups) = [];
    vHe2_tmp(stime_dups) = [];
    vC5_tmp(stime_dups) = [];
    vO6_tmp(stime_dups) = [];
    vFe10_tmp(stime_dups) = [];
    vthHe2_tmp(stime_dups)  = [];
    vthC5_tmp(stime_dups)  = [];
    vthO6_tmp(stime_dups) = [];
    vthFe10_tmp(stime_dups)= [];
    He_qual_tmp(stime_dups) = [];
    C5_qual_tmp(stime_dups) = [];
    O6_qual_tmp(stime_dups) = [];
    Fe10_qual_tmp(stime_dups) = [];
    C6to5_tmp(stime_dups) = [];
    O7to6_tmp(stime_dups) = [];
    avqC_tmp(stime_dups) = [];
    avqO_tmp(stime_dups) = [];
    avqFe_tmp(stime_dups) = [];
    C6to5_qual_tmp(stime_dups) = [];
    O7to6_qual_tmp(stime_dups) = [];
    avqC_qual_tmp(stime_dups) = [];
    avqO_qual_tmp(stime_dups) = [];
    avqFe_qual_tmp(stime_dups) = [];
    SW_type_tmp(stime_dups) = [];
    FetoO_tmp(stime_dups) = [];
    FetoO_qual_tmp(stime_dups) = [];

    stime_skips = find(~(diff(stime_tmp) > (1/24-1/1440) & diff(stime_tmp) < (1/24+1/1440) ) ); # match +/1 1 minute of hour
    numel(stime_skips);
    fflush(stdout)
    for j=1:numel(stime_skips)

      stime_tmp = [stime_tmp(1:stime_skips(j)); stime_tmp(stime_skips(j))+1/24; stime_tmp(stime_skips(j)+1:end)];

      nHe2_tmp = [nHe2_tmp(1:stime_skips(j)); -9999.9; nHe2_tmp(stime_skips(j)+1:end)];
      vHe2_tmp = [vHe2_tmp(1:stime_skips(j)); -9999.9; vHe2_tmp(stime_skips(j)+1:end)];
      vC5_tmp = [vC5_tmp(1:stime_skips(j)); -9999.9; vC5_tmp(stime_skips(j)+1:end)];
      vO6_tmp = [vO6_tmp(1:stime_skips(j)); -9999.9; vO6_tmp(stime_skips(j)+1:end)];
      vFe10_tmp = [vFe10_tmp(1:stime_skips(j)); -9999.9; vFe10_tmp(stime_skips(j)+1:end)];
      vthHe2_tmp = [vthHe2_tmp(1:stime_skips(j)); -9999.9; vthHe2_tmp(stime_skips(j)+1:end)];
      vthC5_tmp = [vthC5_tmp(1:stime_skips(j)); -9999.9; vthC5_tmp(stime_skips(j)+1:end)];
      vthO6_tmp = [vthO6_tmp(1:stime_skips(j)); -9999.9; vthO6_tmp(stime_skips(j)+1:end)];
      vthFe10_tmp = [vthFe10_tmp(1:stime_skips(j)); -9999.9; vthFe10_tmp(stime_skips(j)+1:end)];
      He_qual_tmp = [He_qual_tmp(1:stime_skips(j)); -9999.9; He_qual_tmp(stime_skips(j)+1:end)];
      C5_qual_tmp = [C5_qual_tmp(1:stime_skips(j)); -9999.9; C5_qual_tmp(stime_skips(j)+1:end)];
      O6_qual_tmp = [O6_qual_tmp(1:stime_skips(j)); -9999.9; O6_qual_tmp(stime_skips(j)+1:end)];
      Fe10_qual_tmp = [Fe10_qual_tmp(1:stime_skips(j)); -9999.9; Fe10_qual_tmp(stime_skips(j)+1:end)];
      C6to5_tmp = [C6to5_tmp(1:stime_skips(j)); -9999.9; C6to5_tmp(stime_skips(j)+1:end)];
      O7to6_tmp = [O7to6_tmp(1:stime_skips(j)); -9999.9; O7to6_tmp(stime_skips(j)+1:end)];
      avqC_tmp = [avqC_tmp(1:stime_skips(j)); -9999.9; avqC_tmp(stime_skips(j)+1:end)];
      avqO_tmp = [avqO_tmp(1:stime_skips(j)); -9999.9; avqO_tmp(stime_skips(j)+1:end)];
      avqFe_tmp = [avqFe_tmp(1:stime_skips(j)); -9999.9; avqFe_tmp(stime_skips(j)+1:end)];
      C6to5_qual_tmp = [C6to5_qual_tmp(1:stime_skips(j)); -9999.9; C6to5_qual_tmp(stime_skips(j)+1:end)];
      O7to6_qual_tmp = [O7to6_qual_tmp(1:stime_skips(j)); -9999.9; O7to6_qual_tmp(stime_skips(j)+1:end)];
      avqC_qual_tmp = [avqC_qual_tmp(1:stime_skips(j)); -9999.9; avqC_qual_tmp(stime_skips(j)+1:end)];
      avqO_qual_tmp = [avqO_qual_tmp(1:stime_skips(j)); -9999.9; avqO_qual_tmp(stime_skips(j)+1:end)];
      avqFe_qual_tmp = [avqFe_qual_tmp(1:stime_skips(j)); -9999.9; avqFe_qual_tmp(stime_skips(j)+1:end)];
      SW_type_tmp = [SW_type_tmp(1:stime_skips(j)); -9999.9; SW_type_tmp(stime_skips(j)+1:end)];
      FetoO_tmp = [FetoO_tmp(1:stime_skips(j)); -9999.9; FetoO_tmp(stime_skips(j)+1:end)];
      FetoO_qual_tmp = [FetoO_qual_tmp(1:stime_skips(j)); -9999.9; FetoO_qual_tmp(stime_skips(j)+1:end)];


      stime_skips(j:end) = stime_skips(j:end) + 1;

    endfor


    # append the vectors, then (re)sort by serial time
    [stime_swics, idx] = sort([stime_swics; stime_tmp]);
    nHe2 = [nHe2; nHe2_tmp](idx);
    vHe2 = [vHe2; vHe2_tmp](idx);
    vC5 = [vC5; vC5_tmp](idx);
    vO6 = [vO6; vO6_tmp](idx);
    vFe10 = [vFe10; vFe10_tmp](idx);
    vthHe2 = [vthHe2; vthHe2_tmp](idx);
    vthC5 = [vthC5; vthC5_tmp](idx);
    vthO6 = [vthO6; vthO6_tmp](idx);
    vthFe10 = [vthFe10; vthFe10_tmp](idx);
    He_qual = [He_qual; He_qual_tmp](idx);
    C5_qual = [C5_qual; C5_qual_tmp](idx);
    O6_qual = [O6_qual; O6_qual_tmp](idx);
    Fe10_qual = [Fe10_qual; Fe10_qual_tmp](idx);
    C6to5 = [C6to5; C6to5_tmp](idx);
    O7to6 = [O7to6; O7to6_tmp](idx);
    avqC = [avqC; avqC_tmp](idx);
    avqO = [avqO; avqO_tmp](idx);
    avqFe = [avqFe; avqFe_tmp](idx);
    C6to5_qual = [C6to5_qual; C6to5_qual_tmp](idx);
    O7to6_qual = [O7to6_qual; O7to6_qual_tmp](idx);
    avqC_qual = [avqC_qual; avqC_qual_tmp](idx);
    avqO_qual = [avqO_qual; avqO_qual_tmp](idx);
    avqFe_qual = [avqFe_qual; avqFe_qual_tmp](idx);
    SW_type = [SW_type; SW_type_tmp](idx);
    FetoO = [FetoO; FetoO_tmp];
    FetoO_qual = [FetoO_qual; FetoO_qual_tmp];



    printf("\rDone reading %s (file %d of %d)", swics_files{i}, i, numel(swics_files));
    fflush(stdout);

  endfor


  # convert bad measurements to NaNs for subsequent processing
  nHe2(nHe2 == -9999.9) = nan;
  vHe2(vHe2 == -9999.9) = nan;
  vC5(vC5 == -9999.9) = nan;
  vO6(vO6 == -9999.9) = nan;
  vFe10(vFe10 == -9999.9) = nan;
  vthHe2(vthHe2 == -9999.9) = nan;
  vthC5(vthC5 == -9999.9) = nan;
  vthO6(vthO6 == -9999.9) = nan;
  vthFe10(vthFe10 == -9999.9) = nan;
  He_qual(He_qual == -9999.9) = nan;
  C5_qual(C5_qual == -9999.9) = nan;
  O6_qual(O6_qual == -9999.9) = nan;
  Fe10_qual(Fe10_qual == -9999.9) = nan;
  C6to5(C6to5 == -9999.9) = nan;
  O7to6(O7to6 == -9999.9) = nan;
  avqC(avqC == -9999.9) = nan;
  avqO(avqO == -9999.9) = nan;
  avqFe(avqFe == -9999.9) = nan;
  SW_type(SW_type == -9999.9) = nan;
  FetoO(FetoO == -9999.9) = nan;
  FetoO_qual(FetoO_qual == -9999.9) = nan;



  printf("\n");
  clear *_tmp i j hdr_lines idx fid



endif


#
# Construct multichannel inputs, train classifier, generate output
#
if (0)

  # common date range
  stime_beg = datenum(1999,1,1);
  stime_end = datenum(2008,12,31)
  stime_7chan = stime_swics(stime_swics >= stime_beg & stime_swics <= stime_end);

  # FIXME: the ace_test.mat file contains the variable Htemp_exp, but I don't
  #        remember how it was calculated. I think it is a simple empirical
  #        function of distance from the sun, and other SW parameters, but I
  #        cannot be certain. -EJR

  # 7-channel inputs
  ace_7chan = reshape([Hrho(stime_swe >= stime_beg & stime_swe <= stime_end);
                       Htemp_exp(stime_swe >= stime_beg & stime_swe <= stime_end) ./ Htemp(stime_swe >= stime_beg & stime_swe <= stime_end);
                       HeHratio(stime_swe >= stime_beg & stime_swe <= stime_end);
                       Hspeed(stime_swe >= stime_beg & stime_swe <= stime_end);
                       Bmag(stime_mag >= stime_beg & stime_mag <= stime_end);
                       O7to6(stime_swics >= stime_beg & stime_swics <= stime_end);
                       avqFe(stime_swics >= stime_beg & stime_swics <= stime_end)], 87649, 1, 7);

  # Log10 of 7-channel inputs
  ace_7chan_log10 = log10(ace_7chan);

  # Channel labels
  chan_labels = {'N_p (cm^{-3})', 'T_{exp}/T', '\alpha/p', '|V_p| (km s^{-1})', '|B| (nT)', 'O^{+7}/O^{+6}', 'Q_{Fe}'};

  # train the classifier
  [M_7c_4t, C_7c_4t, N_4t, T_4t] = train_class (ace_7chan, 4);

  # generate thematic time series
  [I_4t, P_4t] = icm_map_class(ace_7chan, [], M_7c_4t, C_7c_4t, N_4t, 0, 1.5, 20);

  # generate a mask for plotting
  I_4t_plotmask = nan([size(I_4t),4]);
  I_4t_plotmask(I_4t==1,:,1) = 1;
  I_4t_plotmask(I_4t==2,:,2) = 1;
  I_4t_plotmask(I_4t==3,:,3) = 1;
  I_4t_plotmask(I_4t==4,:,4) = 1;

  # Theme labels
  t_labels = {'slow SW', 'fast SW', 'CME', 'CIR', 'RR'};


endif


#
# Plot some stuff
#

# Plot fitted marginal distributions for each theme for each channel
if (0)

   %figure('visible','off');
   set(gcf,'paperposition', [.5 .5 7.5 10]);

   for k=1:7

      subplot(7,1,k);

      % define min/max of plot range for kth channel
      xmin = min(M_7c_4t(1,k,:) - 3*sqrt(C_7c_4t(k,k,:)));
      xmax = max(M_7c_4t(1,k,:) + 3*sqrt(C_7c_4t(k,k,:)));

      % define 1001 "independent" variables within min/max plot range
      x = [xmin:(xmax-xmin)/1000:xmax]';

      % plot color-coded marginal normal pdf for each theme for kth channel

      plot(x, [normpdf(x, M_7c_4t(1,k,1), sqrt(C_7c_4t(k,k,1))) * (xmax-xmin)/1000, ...
               normpdf(x, M_7c_4t(1,k,2), sqrt(C_7c_4t(k,k,2))) * (xmax-xmin)/1000, ...
               normpdf(x, M_7c_4t(1,k,3), sqrt(C_7c_4t(k,k,3))) * (xmax-xmin)/1000, ...
               normpdf(x, M_7c_4t(1,k,4), sqrt(C_7c_4t(k,k,4))) * (xmax-xmin)/1000 ], ...
           'linewidth', 4);

      axis([xmin, xmax]);

      %xlabel();
      ylabel(chan_labels{k});
      legend(t_labels);

   endfor % for k=1:7

   % create an image file
   print('-dpdf', 'ACE_7channels_4themes_MargPDFs.pdf');

endif


# Plot 5 1-year panels on a page
if (0)

  # create stime index arrays for years 1999-2008
  idx_1999_7chan = find(stime_7chan >= datenum(1999,1,1) & stime_7chan < datenum(2000,1,1));
  idx_2000_7chan = find(stime_7chan >= datenum(2000,1,1) & stime_7chan < datenum(2001,1,1));
  idx_2001_7chan = find(stime_7chan >= datenum(2001,1,1) & stime_7chan < datenum(2002,1,1));
  idx_2002_7chan = find(stime_7chan >= datenum(2002,1,1) & stime_7chan < datenum(2003,1,1));
  idx_2003_7chan = find(stime_7chan >= datenum(2003,1,1) & stime_7chan < datenum(2004,1,1));
  idx_2004_7chan = find(stime_7chan >= datenum(2004,1,1) & stime_7chan < datenum(2005,1,1));
  idx_2005_7chan = find(stime_7chan >= datenum(2005,1,1) & stime_7chan < datenum(2006,1,1));
  idx_2006_7chan = find(stime_7chan >= datenum(2006,1,1) & stime_7chan < datenum(2009,1,1));
  idx_2006_7chan = find(stime_7chan >= datenum(2006,1,1) & stime_7chan < datenum(2007,1,1));
  idx_2007_7chan = find(stime_7chan >= datenum(2007,1,1) & stime_7chan < datenum(2008,1,1));
  idx_2008_7chan = find(stime_7chan >= datenum(2008,1,1) & stime_7chan < datenum(2009,1,1));

  # plot each year as 1 of 5 panels per page
  figure(1);
  set(gcf, 'paperposition', [.5 .5 7.5 10]);

  subplot(5,1,1);
  plot(stime_7chan(idx_1999_7chan), repmat(ace_7chan(idx_1999_7chan,:,4),1,1), 'k');
  axis([datenum(1999,1,1) datenum(1999,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_1999_7chan), repmat(ace_7chan(idx_1999_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_1999_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,2);
  plot(stime_7chan(idx_2000_7chan), repmat(ace_7chan(idx_2000_7chan,:,4),1,1), 'k');
  axis([datenum(2000,1,1) datenum(2000,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2000_7chan), repmat(ace_7chan(idx_2000_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2000_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,3);
  plot(stime_7chan(idx_2001_7chan), repmat(ace_7chan(idx_2001_7chan,:,4),1,1), 'k');
  axis([datenum(2001,1,1) datenum(2001,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2001_7chan), repmat(ace_7chan(idx_2001_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2001_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,4);
  plot(stime_7chan(idx_2002_7chan), repmat(ace_7chan(idx_2002_7chan,:,4),1,1), 'k');
  axis([datenum(2002,1,1) datenum(2002,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2002_7chan), repmat(ace_7chan(idx_2002_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2002_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,5);
  plot(stime_7chan(idx_2003_7chan), repmat(ace_7chan(idx_2003_7chan,:,4),1,1), 'k');
  axis([datenum(2003,1,1) datenum(2003,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2003_7chan), repmat(ace_7chan(idx_2003_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2003_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  print('-dpdf','Vsw_4themes_1999to2003.pdf')


  figure(2);
  set(gcf, 'paperposition', [.5 .5 7.5 10]);

  subplot(5,1,1);
  plot(stime_7chan(idx_2004_7chan), repmat(ace_7chan(idx_2004_7chan,:,4),1,1), 'k');
  axis([datenum(2004,1,1) datenum(2004,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2004_7chan), repmat(ace_7chan(idx_2004_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2004_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,2);
  plot(stime_7chan(idx_2005_7chan), repmat(ace_7chan(idx_2005_7chan,:,4),1,1), 'k');
  axis([datenum(2005,1,1) datenum(2005,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2005_7chan), repmat(ace_7chan(idx_2005_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2005_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,3);
  plot(stime_7chan(idx_2006_7chan), repmat(ace_7chan(idx_2006_7chan,:,4),1,1), 'k');
  axis([datenum(2006,1,1) datenum(2006,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2006_7chan), repmat(ace_7chan(idx_2006_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2006_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,4);
  plot(stime_7chan(idx_2007_7chan), repmat(ace_7chan(idx_2007_7chan,:,4),1,1), 'k');
  axis([datenum(2007,1,1) datenum(2007,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2007_7chan), repmat(ace_7chan(idx_2007_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2007_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  subplot(5,1,5);
  plot(stime_7chan(idx_2008_7chan), repmat(ace_7chan(idx_2008_7chan,:,4),1,1), 'k');
  axis([datenum(2008,1,1) datenum(2008,12,31) 200 1200]);
  hold on;
  plot(stime_7chan(idx_2008_7chan), repmat(ace_7chan(idx_2008_7chan,:,4),1,4) .* squeeze(I_4t_plotmask(idx_2008_7chan,:,:)), 'linewidth', 4);
  datetick('x', 'yyyy-mm-dd', 'keeplimits');
  grid on;
  legend('','slow SW','fast SW','CME','CIR');
  ylabel('km/s');

  print('-dpdf','Vsw_4themes_2004to2008.pdf')


endif


# Plot 4 3-month panels on a page
if (1)

  idx_1999q1_7chan = find(stime_7chan >= datenum(1999,1,1) & stime_7chan < datenum(1999,4,1));
  idx_1999q2_7chan = find(stime_7chan >= datenum(1999,4,1) & stime_7chan < datenum(1999,7,1));
  idx_1999q3_7chan = find(stime_7chan >= datenum(1999,7,1) & stime_7chan < datenum(1999,10,1));
  idx_1999q4_7chan = find(stime_7chan >= datenum(1999,10,1) & stime_7chan < datenum(2000,1,1));

  idx_2000q1_7chan = find(stime_7chan >= datenum(2000,1,1) & stime_7chan < datenum(2000,4,1));
  idx_2000q2_7chan = find(stime_7chan >= datenum(2000,4,1) & stime_7chan < datenum(2000,7,1));
  idx_2000q3_7chan = find(stime_7chan >= datenum(2000,7,1) & stime_7chan < datenum(2000,10,1));
  idx_2000q4_7chan = find(stime_7chan >= datenum(2000,10,1) & stime_7chan < datenum(2001,1,1));

  idx_2001q1_7chan = find(stime_7chan >= datenum(2001,1,1) & stime_7chan < datenum(2001,4,1));
  idx_2001q2_7chan = find(stime_7chan >= datenum(2001,4,1) & stime_7chan < datenum(2001,7,1));
  idx_2001q3_7chan = find(stime_7chan >= datenum(2001,7,1) & stime_7chan < datenum(2001,10,1));
  idx_2001q4_7chan = find(stime_7chan >= datenum(2001,10,1) & stime_7chan < datenum(2002,1,1));

  idx_2002q1_7chan = find(stime_7chan >= datenum(2002,1,1) & stime_7chan < datenum(2002,4,1));
  idx_2002q2_7chan = find(stime_7chan >= datenum(2002,4,1) & stime_7chan < datenum(2002,7,1));
  idx_2002q3_7chan = find(stime_7chan >= datenum(2002,7,1) & stime_7chan < datenum(2002,10,1));
  idx_2002q4_7chan = find(stime_7chan >= datenum(2002,10,1) & stime_7chan < datenum(2003,1,1));

  idx_2003q1_7chan = find(stime_7chan >= datenum(2003,1,1) & stime_7chan < datenum(2003,4,1));
  idx_2003q2_7chan = find(stime_7chan >= datenum(2003,4,1) & stime_7chan < datenum(2003,7,1));
  idx_2003q3_7chan = find(stime_7chan >= datenum(2003,7,1) & stime_7chan < datenum(2003,10,1));
  idx_2003q4_7chan = find(stime_7chan >= datenum(2003,10,1) & stime_7chan < datenum(2004,1,1));

  idx_2004q1_7chan = find(stime_7chan >= datenum(2004,1,1) & stime_7chan < datenum(2004,4,1));
  idx_2004q2_7chan = find(stime_7chan >= datenum(2004,4,1) & stime_7chan < datenum(2004,7,1));
  idx_2004q3_7chan = find(stime_7chan >= datenum(2004,7,1) & stime_7chan < datenum(2004,10,1));
  idx_2004q4_7chan = find(stime_7chan >= datenum(2004,10,1) & stime_7chan < datenum(2005,1,1));

  idx_2005q1_7chan = find(stime_7chan >= datenum(2005,1,1) & stime_7chan < datenum(2005,4,1));
  idx_2005q2_7chan = find(stime_7chan >= datenum(2005,4,1) & stime_7chan < datenum(2005,7,1));
  idx_2005q3_7chan = find(stime_7chan >= datenum(2005,7,1) & stime_7chan < datenum(2005,10,1));
  idx_2005q4_7chan = find(stime_7chan >= datenum(2005,10,1) & stime_7chan < datenum(2006,1,1));

  idx_2006q1_7chan = find(stime_7chan >= datenum(2006,1,1) & stime_7chan < datenum(2006,4,1));
  idx_2006q2_7chan = find(stime_7chan >= datenum(2006,4,1) & stime_7chan < datenum(2006,7,1));
  idx_2006q3_7chan = find(stime_7chan >= datenum(2006,7,1) & stime_7chan < datenum(2006,10,1));
  idx_2006q4_7chan = find(stime_7chan >= datenum(2006,10,1) & stime_7chan < datenum(2007,1,1));

  idx_2007q1_7chan = find(stime_7chan >= datenum(2007,1,1) & stime_7chan < datenum(2007,4,1));
  idx_2007q2_7chan = find(stime_7chan >= datenum(2007,4,1) & stime_7chan < datenum(2007,7,1));
  idx_2007q3_7chan = find(stime_7chan >= datenum(2007,7,1) & stime_7chan < datenum(2007,10,1));
  idx_2007q4_7chan = find(stime_7chan >= datenum(2007,10,1) & stime_7chan < datenum(2008,1,1));

  idx_2008q1_7chan = find(stime_7chan >= datenum(2008,1,1) & stime_7chan < datenum(2008,4,1));
  idx_2008q2_7chan = find(stime_7chan >= datenum(2008,4,1) & stime_7chan < datenum(2008,7,1));
  idx_2008q3_7chan = find(stime_7chan >= datenum(2008,7,1) & stime_7chan < datenum(2008,10,1));
  idx_2008q4_7chan = find(stime_7chan >= datenum(2008,10,1) & stime_7chan < datenum(2009,1,1));


  # alternative plot formats
  if (1)

    for i=1:10
      for j=1:4

        #figure((i-1)*4 + j);
        figure('visible','off');
        set(gcf,'paperposition', [.5 .5 7.5 10]);

        eval(sprintf('dateticks = datenum(%4d, [1:12]'', 1);', 1998+i));

        for k=1:7

          subplot(7,1,k);

          eval(sprintf('plot(stime_7chan(idx_%4dq%d_7chan), ace_7chan(idx_%4dq%d_7chan,:,%d), ''k'')', 1998+i, j, 1998+i, j, k) );
          hold on;
          eval([sprintf('plot([0;stime_7chan(idx_%4dq%d_7chan)], ', 1998+i, j), ...
                sprintf('[0 0 0 0; repmat(ace_7chan(idx_%4dq%d_7chan,:,%d),1,4)] .* ', 1998+i, j, k), ...
                sprintf('[0 0 0 0; squeeze(I_4t_plotmask(idx_%4dq%d_7chan,:,:))], ''linewidth'', 4)', 1998+i, j) ] );
          yax = axis()(3:4);
          eval(sprintf('axis([floor(stime_7chan(idx_%4dq%d_7chan(1))) ceil(stime_7chan(idx_%4dq%d_7chan(end))) yax(1) yax(2)])', 1998+i, j, 1998+i, j) );
          set(gca,'xtick', dateticks);
          datetick('x','yyyy-mm-dd', 'keeplimits', 'keepticks');
          grid on;
          legend({'', t_labels{:}});
          ylabel(chan_labels{k});

          hold off;

        endfor % k=1:7 ...channels

        #eval(sprintf('print(''-dpdf'',''ACE_7channels_4themes_%4d_Q%d.pdf'')', 1998+i, j) );
        eval(sprintf('print(''-dpdf'',''ACE_7channels_4themes_%4d_Q%d.pdf'')', 1998+i, j) );

      endfor % j=1:4 ...quarters
    endfor % i=1:10 ...years




  else

    # This format looks horrible, but there might be some useful hints below, so leave it in for now

    figure(1);
    set(gcf, 'paperposition', [.5 .5 7.5 10]);

    subplot(1,1,1);
    # plot Vsw time series
    plot(stime_7chan(idx_1999q1_7chan), ace_7chan(idx_1999q1_7chan,:,4), 'k');
    axis([datenum(1999,1,1) datenum(1999,4,1) 200 1200]);
    hold on;

    xax = axis()(1:2);
    yax = axis()(3:4);


    # plot patches for each theme
    for j=1:4

      idx_tmp = find(~(isnan(I_4t_plotmask(idx_1999q1_7chan,:,j))));

      % don't draw a patch for each index, but rather ranges of indices
      first_first = min(min(0,idx_tmp)); % this is a total kludge to ensure an empty matrix or 0
      last_last = min(min(0,idx_tmp)) + length(idx_tmp); % this is a continuation of total kludge
      first_idx = [first_first; find(diff(idx_tmp) ~= 1)] + 1;
      last_idx = [find(diff(idx_tmp) ~= 1); last_last];


      for k=1:numel(first_idx)
        patch([stime_7chan(idx_1999q1_7chan)(idx_tmp(first_idx(k)))-.5/24 ...
               stime_7chan(idx_1999q1_7chan)(idx_tmp(last_idx(k)))+.5/24 ...
               stime_7chan(idx_1999q1_7chan)(idx_tmp(last_idx(k)))+.5/24 ...
               stime_7chan(idx_1999q1_7chan)(idx_tmp(first_idx(k)))-.5/24], ...
              [yax(1)+.01*diff(yax) yax(1)+.01*diff(yax) yax(2)-.01*diff(yax) yax(2)-.01*diff(yax)], ...
              colormap()(j*(rows(colormap)/4),:), 'linewidth', 0);

        idx_tmp_patch = find(stime_7chan(idx_1999q1_7chan) >= stime_7chan(idx_1999q1_7chan)(idx_tmp(first_idx(k))) & ...
                             stime_7chan(idx_1999q1_7chan) <= stime_7chan(idx_1999q1_7chan)(idx_tmp(last_idx(k))) );

        plot(stime_7chan(idx_1999q1_7chan)(idx_tmp_patch), ace_7chan(idx_1999q1_7chan,:,4)(idx_tmp_patch), ...
             'linewidth', 2, 'color', [1 - (colormap()(j*(rows(colormap)/4),:))]*0+1 );

      endfor


    endfor

    datetick('x', 'yyyy-mm-dd', 'keeplimits');
    grid on;
    legend('','slow SW','fast SW','CME','CIR');
    ylabel('km/s');

    hold off


  endif

endif
