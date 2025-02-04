%combine_sites.m
% combine output from 

clear
close all

import('utilities.readDateTime');

fs = 256000;

%Load validated pings

PING.COC = load("D:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\COC_2020_09\COC_EK60_VALIDATED_RANGES.mat");
PING.GBK = load("D:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\GBK_2020_09\GBK_EK60_VALIDATED_RANGES.mat");
PING.FCH = load("D:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\FCH_2020_09\FCH_EK60_VALIDATED_RANGES.mat");

PING_COC = PING.COC.Validated_ranges;
PING_GBK = PING.GBK.Validated_ranges;
PING_FCH = PING.FCH.Validated_ranges;

PING_FCH = removevars(PING_FCH, "VarName1");

PING_COC.SITE = repmat('COC',[height(PING_COC),1]);
PING_GBK.SITE = repmat('GBK',[height(PING_GBK),1]);
PING_FCH.SITE = repmat('FCH',[height(PING_FCH),1]);

PING_COC.horRange = PING_COC.COCranges;
PING_GBK.horRange = PING_GBK.GBKranges;
PING_FCH.horRange = PING_FCH.FCHranges;

PING_COC.slantRange = PING_COC.COC_SLANT_RANGE;
PING_GBK.slantRange = PING_GBK.GBK_SLANT_RANGE;
PING_FCH.slantRange = PING_FCH.FCH_SLANT_RANGE;

ALL_SITES_PINGS = [PING_COC;PING_GBK;PING_FCH];

%simplify data

ALL_SITES_PINGS = removevars(ALL_SITES_PINGS,{'COC_SLANT_RANGE','FCH_SLANT_RANGE','GBK_SLANT_RANGE','COCranges','FCHranges','GBKranges'});

save("D:\BW_ECHO_EXPERIMENT\DETECTOR_OUTPUT\ALL_SITES\ALL_SITES_PINGS.mat",'ALL_SITES_PINGS')






