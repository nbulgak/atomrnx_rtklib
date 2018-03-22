#include <stdio.h>
#include <stdarg.h>
#include "rtklib.h"

double RestorePValue(double Reference, double Ambiguity, double Modulo)
{

	Modulo = Modulo - (int)(Modulo / Ambiguity)*Ambiguity;
	double Closest = (int)(Reference / Ambiguity)*Ambiguity; /*find nearest on ambiguity scale*/
	double Full = Closest + Modulo; /*full value which can have error +-fAmbiguity*/
	double Delta = (Full - Reference); /*Find how far restored value from reference*/
	if (Delta > Ambiguity / 2) Full -= Ambiguity; /*make minus correction*/
	else if (Delta < -Ambiguity / 2) Full += Ambiguity; /*make plus correction*/
	return Full;
}

static int obsindex(obs_t *obs, gtime_t time, int sat)
{
	int i, j;

	for (i = 0; i<obs->n; i++) {
		if (obs->data[i].sat == sat) return i; /* field already exists */
	}
	if (i >= MAXOBS) return -1; /* overflow */

	/* add new field */
	obs->data[i].time = time;
	obs->data[i].sat = sat;
	for (j = 0; j<NFREQ + NEXOBS; j++) {
		obs->data[i].L[j] = obs->data[i].P[j] = 0.0;
		obs->data[i].D[j] = 0.0;
		obs->data[i].SNR[j] = obs->data[i].LLI[j] = obs->data[i].code[j] = 0;
	}
	obs->n++;
	return i;
}

static void sigindex(int sys, const unsigned char *code, const int *freq, int n,
	const char *opt, int *ind)
{
	int i, nex, pri, pri_h[8] = { 0 }, index[8] = { 0 }, ex[32] = { 0 }, j;

	/* test code priority */

	for (i = 0; i<n; i++) {
		if (!code[i]) continue;

		if (freq[i]>NFREQ) { /* save as extended signal if freq > NFREQ */
			ex[i] = 1;
			continue;
		}
		/* code priority */
		pri = getcodepri(sys, code[i], opt);

		/* select highest priority signal */
		if (pri>pri_h[freq[i] - 1]) {
			if (index[freq[i] - 1]) ex[index[freq[i] - 1] - 1] = 1;
			pri_h[freq[i] - 1] = pri;
			index[freq[i] - 1] = i + 1;
		}
		else ex[i] = 1;
	}

	/* signal index in obs data */
	for (i = nex = 0; i<n; i++) {
		if (ex[i] == 0) ind[i] = freq[i] - 1;
		else if (nex<NEXOBS) ind[i] = NFREQ + nex++;
		else { /* no space in obs data */
			trace(2, "rtcm msm: no space in obs data sys=%d code=%d\n", sys, code[i]);
			ind[i] = -1;
		}
#if 0
		trace(2, "sig pos: sys=%d code=%d ex=%d ind=%d\n", sys, code[i], ex[i], ind[i]);
#endif
	}
}

#if 1
static FILE *nik_log;

static int nik_printf(const char *format, ...)
{
	char buf[256];
	void *tbl;
	va_list args;
	unsigned int size;

	va_start(args, format);
	size = vsnprintf(buf, sizeof(buf), format, args);
	va_end(args);

	if (!nik_log)
		nik_log = fopen("nik_log.txt", "wt");

	if (nik_log)
	{
		fprintf(nik_log, "%s", buf);
		fflush(nik_log);
	}

	return 0;
}
#else
static int nik_printf(const char *format, ...)
{
	return 0;
}
#endif
#define MAX_SATS 64
#define MAX_SIGS 64
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define RANGE_MS    (CLIGHT*0.001)      /* range in 1 ms */
#define printf nik_printf

static void outnavf(double value)
{
	double e = fabs(value)<1E-99 ? 0.0 : floor(log10(fabs(value)) + 1.0);
	printf(" %s.%012.0fE%+03.0f", value<0.0 ? "-" : " ", fabs(value) / pow(10.0, e - 12.0), e);
}

static int get_active_bits(unsigned int mask, int len)
{
	int res;
	for (res = 0; mask > 0 && len > 0; mask >>= 1, len--)
	if ((mask & 1))
		res++;
	return(res);
}

static double getbitg(const unsigned char *buff, int pos, int len)
{
	double value = getbitu(buff, pos + 1, len - 1);
	return getbitu(buff, pos, 1) ? -value : value;
}

/* decode type 1019: gps ephemerides -----------------------------------------*/
static int decode_type1019(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	eph_t eph = { 0 };
	double toc, sqrtA;
	char *msg;
	int i = k, prn, sat, week, sys = SYS_GPS, www;
	double sec_in_week;
	printf("offset=%d\n", k);

	if (i + 476 <= mes_len * 8) {
		prn = getbitu(Raw, i, 6);              i += 6;
		printf("prn=%d\n", prn);
		week = getbitu(Raw, i, 10);              i += 10;
		printf("week=%d\n", week);
		eph.sva = getbitu(Raw, i, 4);              i += 4;
		printf("sva=%d\n", eph.sva);
		eph.code = getbitu(Raw, i, 2);              i += 2;
		printf("code=%d\n", eph.code);
		eph.idot = getbits(Raw, i, 14)*P2_43*SC2RAD; i += 14;
		printf("idot=%.15lf\n", eph.idot);
		eph.iode = getbitu(Raw, i, 8);              i += 8;
		printf("iode=%d\n", eph.iode);
		toc = getbitu(Raw, i, 16)*16.0;         i += 16;
		printf("toc=%.15lf\n", toc);
		eph.f2 = getbits(Raw, i, 8)*P2_55;        i += 8;
		printf("f2=%.15lf\n", eph.f2);
		eph.f1 = getbits(Raw, i, 16)*P2_43;        i += 16;
		printf("f1=%.15lf\n", eph.f1);
		eph.f0 = getbits(Raw, i, 22)*P2_31;        i += 22;
		printf("f0=%.15lf\n", eph.f0);
		eph.iodc = getbitu(Raw, i, 10);              i += 10;
		printf("iodc=%d\n", eph.iodc);
		eph.crs = getbits(Raw, i, 16)*P2_5;         i += 16;
		printf("crs=%.15lf\n", eph.crs);
		eph.deln = getbits(Raw, i, 16)*P2_43*SC2RAD; i += 16;
		printf("deln=%.15lf\n", eph.deln);
		eph.M0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		printf("M0=%.15lf\n", eph.M0);
		eph.cuc = getbits(Raw, i, 16)*P2_29;        i += 16;
		printf("cuc=%.15lf\n", eph.cuc);
		eph.e = getbitu(Raw, i, 32)*P2_33;        i += 32;
		printf("e=%.15lf\n", eph.e);
		eph.cus = getbits(Raw, i, 16)*P2_29;        i += 16;
		printf("cus=%.15lf\n", eph.cus);
		sqrtA = getbitu(Raw, i, 32)*P2_19;        i += 32;
		printf("sqrtA=%.15lf\n", sqrtA);
		eph.toes = getbitu(Raw, i, 16)*16.0;         i += 16;
		printf("toes=%.15lf\n", eph.toes);
		eph.cic = getbits(Raw, i, 16)*P2_29;        i += 16;
		printf("cic=%.15lf\n", eph.cic);
		eph.OMG0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		printf("OMG0=%.15lf\n", eph.OMG0);
		eph.cis = getbits(Raw, i, 16)*P2_29;        i += 16;
		printf("cis=%.15lf\n", eph.cis);
		eph.i0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		printf("i0=%.15lf\n", eph.i0);
		eph.crc = getbits(Raw, i, 16)*P2_5;         i += 16;
		printf("crc=%.15lf\n", eph.crc);
		eph.omg = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		printf("omg=%.15lf\n", eph.omg);
		eph.OMGd = getbits(Raw, i, 24)*P2_43*SC2RAD; i += 24;
		printf("OMGd=%.15lf\n", eph.OMGd);
		eph.tgd[0] = getbits(Raw, i, 8)*P2_31;        i += 8;
		printf("tgd[0]=%.15lf\n", eph.tgd[0]);
		eph.svh = getbitu(Raw, i, 6);              i += 6;
		printf("svh=%d\n", eph.svh);
		eph.flag = getbitu(Raw, i, 1);              i += 1;
		printf("flag=%d\n", eph.flag);
		eph.fit = getbitu(Raw, i, 1) ? 0.0 : 4.0; /* 0:4hr,1:>4hr */
		printf("fit=%d\n", eph.fit);
	}
	else {
		printf("Invalid message: Len_NAV_1019\n");
		return -1;
	}
	if (prn >= 40) {
		/*sys = SYS_SBS; prn += 80;*/
		return 0; /*udalit kogda nujen budet SBS*/
	}
	/*	trace(4, "decode_type1019: prn=%d iode=%d toe=%.0f\n", prn, eph.iode, eph.toes);

		if (rtcm->outtype) {
		msg = rtcm->msgtype + strlen(rtcm->msgtype);
		sprintf(msg, " prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
		prn, eph.iode, eph.iodc, week, eph.toes, toc, eph.svh);
		}
		*/
	if (!(sat = satno(sys, prn))) {
		return -1;
	}
	eph.sat = sat;
	printf("SAT=%d\n", eph.sat);
	eph.week = adjgpsweek(week);

    sec_in_week = time2gpst(raw->time, &www);
	raw->time = gpst2time(eph.week, sec_in_week);

	printf("WEEK=%d\n", eph.week);
	eph.toe = gpst2time(eph.week, eph.toes);
	eph.toc = gpst2time(eph.week, toc);
	eph.A = sqrtA*sqrtA;

   	if (eph.iode == raw->nav.eph[sat - 1].iode) return 0; /* unchanged */


	/*printf("number of broadcast ephemeris. n=%d, nmax=%d\n", raw->nav.n, raw->nav.nmax);*/

	raw->nav.eph[sat - 1] = eph;
	raw->ephsat = sat;
	printf("SAT[%d] dobavlen\n", sat);
	return 2;
}
/* decode type 1020: glonass ephemerides -------------------------------------*/
static int decode_type1020(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	geph_t geph = { 0 };
	double tk_h, tk_m, tk_s, toe, tow, tod, tof;
	char *msg, s[64];
	int i = k, prn, sat, week, tb, bn, sys = SYS_GLO;

	if (i + 348 <= mes_len * 8) {
		prn = getbitu(Raw, i, 6);           i += 6;
		outnavf(prn); printf("\n");
		geph.frq = getbitu(Raw, i, 5) - 7;         i += 5 + 2 + 2;
		outnavf(geph.frq); printf("\n");
		tk_h = getbitu(Raw, i, 5);           i += 5;
		outnavf(tk_h); printf("\n");
		tk_m = getbitu(Raw, i, 6);           i += 6;
		outnavf(tk_m); printf("\n");
		tk_s = getbitu(Raw, i, 1)*30.0;      i += 1;
		outnavf(tk_s); printf("\n");
		bn = getbitu(Raw, i, 1);           i += 1 + 1;
		outnavf(bn); printf("\n");
		tb = getbitu(Raw, i, 7);           i += 7;
		outnavf(tb); printf("\n");
		geph.vel[0] = getbitg(Raw, i, 24)*P2_20*1E3; i += 24;
		outnavf(geph.vel[0]); printf("\n");
		geph.pos[0] = getbitg(Raw, i, 27)*P2_11*1E3; i += 27;
		outnavf(geph.pos[0]); printf("\n");
		geph.acc[0] = getbitg(Raw, i, 5)*P2_30*1E3; i += 5;
		outnavf(geph.acc[0]); printf("\n");
		geph.vel[1] = getbitg(Raw, i, 24)*P2_20*1E3; i += 24;
		outnavf(geph.vel[1]); printf("\n");
		geph.pos[1] = getbitg(Raw, i, 27)*P2_11*1E3; i += 27;
		outnavf(geph.pos[1]); printf("\n");
		geph.acc[1] = getbitg(Raw, i, 5)*P2_30*1E3; i += 5;
		outnavf(geph.acc[1]); printf("\n");
		geph.vel[2] = getbitg(Raw, i, 24)*P2_20*1E3; i += 24;
		outnavf(geph.vel[2]); printf("\n");
		geph.pos[2] = getbitg(Raw, i, 27)*P2_11*1E3; i += 27;
		outnavf(geph.pos[2]); printf("\n");
		geph.acc[2] = getbitg(Raw, i, 5)*P2_30*1E3; i += 5 + 1;
		outnavf(geph.acc[2]); printf("\n");
		geph.gamn = getbitg(Raw, i, 11)*P2_40;     i += 11 + 3;
		outnavf(geph.gamn); printf("\n");
		geph.taun = getbitg(Raw, i, 22)*P2_30;
		outnavf(geph.taun); printf("\n");
	}
	else {
		printf("Invalid message: Len_NAV_1020\n");
		/*trace(2, "rtcm3 1020 length error: len=%d\n", rtcm->len);*/
		return -1;
	}
	if (!(sat = satno(sys, prn))) {
		/*trace(2, "rtcm3 1020 satellite number error: prn=%d\n", prn);*/
		return -1;
	}
	/*trace(4, "decode_type1020: prn=%d tk=%02.0f:%02.0f:%02.0f\n", prn, tk_h, tk_m, tk_s);*/

	/*if (rtcm->outtype) {
		msg = rtcm->msgtype + strlen(rtcm->msgtype);
		sprintf(msg, " prn=%2d tk=%02.0f:%02.0f:%02.0f frq=%2d bn=%d tb=%d",
			prn, tk_h, tk_m, tk_s, geph.frq, bn, tb);
	}*/
	geph.sat = sat;
	geph.svh = bn;
	geph.iode = tb & 0x7F;
	/*if (raw->time.time == 0) { printf("ya v raw->time.time == 0\n"); raw->time = utc2gpst(timeget()); }*/
	tow = time2gpst(raw->time, &week);
	printf("week=%d\n", week);
	tod = fmod(tow, 86400.0); tow -= tod;
	tof = tk_h*3600.0 + tk_m*60.0 + tk_s - 10800.0; /* lt->utc */
	if (tof<tod - 43200.0) tof += 86400.0;
	else if (tof>tod + 43200.0) tof -= 86400.0;
	geph.tof = utc2gpst(gpst2time(week, tow + tof));
	toe = tb*900.0 - 10800.0; /* lt->utc */
	if (toe<tod - 43200.0) toe += 86400.0;
	else if (toe>tod + 43200.0) toe -= 86400.0;
	geph.toe = utc2gpst(gpst2time(week, tow + toe)); /* utc->gpst */

	/*if (!strstr(rtcm->opt, "-EPHALL")) {*/
		if (fabs(timediff(geph.toe, raw->nav.geph[prn - 1].toe))<1.0&&
			geph.svh == raw->nav.geph[prn - 1].svh) return 0; /* unchanged */
	/*}*/
	raw->nav.geph[prn - 1] = geph;
	raw->ephsat = sat;
	printf("SAT[%d] dobavlen\n", sat);
	return 2;
}

static int decode_atom_rnx(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	extern const char *msm_sig_gps[32], *msm_sig_glo[32];
	char www[64];
	int i, j, l, m, week;
	double FinePseudoPhase;
	int Doppler[MAX_SATS], FineDoppler[MAX_SATS][MAX_SIGS];
	unsigned int GNSSmask, version, ref_st_ID, multiple, IODS, smoothing, position_presentation, primary_GNSS_system,
		primary_time_tag, time_tag_extension_type, fractional_second, hour, day, divergence_free_smoothing_indicator,
		cumulative_session_transmitting_time_indicator, Data_ID_change_counter, Data_ID_follow, Nms_follow, Supplementary_follow,
		Pseudo_range_follow, Carrier_phase_follow, Resolution, Reserved1, Reserved2, Satellite_mask[MAX_SATS / 32], Signal_mask,
		Satellite_mask_len[MAX_SATS / 32], Signal_mask_len, sat_cnt, sig_cnt, sat[MAX_SATS], sat_sig_cnt[MAX_SATS], sig[MAX_SIGS],
		sig_type[MAX_SATS][MAX_SIGS], offs, ncell, need_bits, Int_num_sat_ranges[MAX_SATS], Sat_rough_range[MAX_SATS], Azimuth[MAX_SATS],
		Elevation[MAX_SATS], Full_Range_Available[MAX_SATS], Satellite_Usage_Status[MAX_SATS], ChannelNumber[MAX_SATS][MAX_SIGS],
		FinePseudoRange[MAX_SATS][MAX_SIGS], CycSlipCounter[MAX_SATS][MAX_SIGS], IntCycPhase[MAX_SATS][MAX_SIGS],
		FracCycPhase[MAX_SATS][MAX_SIGS], SNR[MAX_SATS][MAX_SIGS], ExtSuppData[MAX_SATS][MAX_SIGS][2], Reference_P[MAX_SATS];

	if (k + 58 > mes_len * 8)
	{
		printf("Invalid message: Len_Message Header \n");
		return -1;
	}

	printf("===============================================================RNX\n");
	/*MESSAGE HEADER*/
	version = getbitu(Raw, k, 3); k += 3;
	printf("version=%d\n", version);

	ref_st_ID = getbitu(Raw, k, 12); k += 12;
	printf("reference stantion ID=%d\n", ref_st_ID);

	multiple = getbitu(Raw, k, 1); k += 1;
	printf("multiple message bit=%d\n", multiple);

	IODS = getbitu(Raw, k, 3); k += 3;
	printf("IODS=%d\n", IODS);

	smoothing = getbitu(Raw, k, 3); k += 3;
	printf("Smoothing interval=%d\n", smoothing);

	position_presentation = getbitu(Raw, k, 2); k += 2;
	printf("Position presentation=%d\n", position_presentation);

	GNSSmask = getbitu(Raw, k, 8); k += 8;
	printf("GNSS mask=");
	for (j = 1; j < 129;j*=2)
	    printf("%d", !(!(GNSSmask&j)) );
	printf("\n");
	if (0 == GNSSmask)
	{
		printf("no info message\n");
		return 1;
	}

	primary_GNSS_system = getbitu(Raw, k, 3); k += 3;
	printf("Primary GNSS system=%d\n", primary_GNSS_system);

	primary_time_tag = getbitu(Raw, k, 12); k += 12;
	printf("Primary time tag=%d\n", primary_time_tag);

	time_tag_extension_type = getbitu(Raw, k, 1); k += 1;
	printf("Time tag extension type=%d\n", time_tag_extension_type);

	if (time_tag_extension_type == 1)
	{
		fractional_second = getbitu(Raw, k, 8); k += 8;
		printf("Fractional second=%d\n", fractional_second);
	}
	else
	{
		hour = getbitu(Raw, k, 5); k += 5;
		printf("Hour=%d\n", hour);
		day = getbitu(Raw, k, 3); k += 3;
		printf("Day=%d\n", day);
	}

	divergence_free_smoothing_indicator = getbitu(Raw, k, 1); k++;
	printf("Divergence free smoothing indicator=%d\n", divergence_free_smoothing_indicator);

	cumulative_session_transmitting_time_indicator = getbitu(Raw, k, 7); k += 7;
	printf("Cumulative session transmitting time indicator=%d\n", cumulative_session_transmitting_time_indicator);

	/*Observables mask*/
	if (k + 16 > mes_len * 8)
	{
		printf("Invalid message: Len_Observable mask\n");
		return -1;
	}

	Data_ID_change_counter = getbitu(Raw, k, 5); k += 5;
	printf("Data_ID_change_counter=%d\n", Data_ID_change_counter);

	Data_ID_follow = getbitu(Raw, k, 1); k += 1;
	printf("Data_ID_follow=%d\n", Data_ID_follow);

	Nms_follow = getbitu(Raw, k, 1); k += 1;
	printf("Nms_follow=%d\n", Nms_follow);

	Supplementary_follow = getbitu(Raw, k, 2); k += 2;
	printf("Supplementary_follow=%d\n", Supplementary_follow);

	Pseudo_range_follow = getbitu(Raw, k, 2); k += 2;
	printf("Pseudo_range_follow=%d\n", Pseudo_range_follow);

	Carrier_phase_follow = getbitu(Raw, k, 2); k += 2;
	printf("Carrier_phase_follow=%d\n", Carrier_phase_follow);

	Resolution = getbitu(Raw, k, 1); k += 1;
	printf("Resolution=%d\n", Resolution);

	Reserved1 = getbitu(Raw, k, 1); k += 1;
	printf("Reserved1=%d\n", Reserved1);

	Reserved2 = getbitu(Raw, k, 1); k += 1;
	printf("Reserved2=%d\n", Reserved2);

	/*Sat/Sig mask*/
	if (Data_ID_follow != 1)
	{
		printf("Invalid message: Data Id follow\n");
		return -1;
	}
	if (version == 2)
	{
		if (k + 96 > mes_len * 8)
		{
			printf("Invalid message: Len_Sat/Sig MASK1\n");
			return -1;
		}

		Satellite_mask[0] = getbitu(Raw, k, 32); k += 32;
		Satellite_mask[1] = getbitu(Raw, k, 32); k += 32;
		printf("Satellite_mask=%08x%08x\n", Satellite_mask[0], Satellite_mask[1]);
		Satellite_mask_len[0] = 32; Satellite_mask_len[1] = 32;

		Signal_mask = getbitu(Raw, k, 32); k += 32;
		printf("Signal_mask=%x\n", Signal_mask);
		Signal_mask_len = 32;
	}
	else
	{
		if (k + 72 > mes_len * 8)
		{
			printf("Invalid message: Len_Sat/Sig MASK2\n");
			return -1;
		}

		Satellite_mask[0] = getbitu(Raw, k, 32); k += 32;
		Satellite_mask[1] = getbitu(Raw, k, 8); k += 8;
		printf("Satellite_mask=%08x%02x\n", Satellite_mask[0], Satellite_mask[1]);
		Satellite_mask_len[0] = 32; Satellite_mask_len[1] = 8;

		Signal_mask = getbitu(Raw, k, 24); k += 24;
		printf("Signal_mask=%x\n", Signal_mask);
		Signal_mask_len = 24;
		k += 8; /* reserved*/
	}

	sat_cnt = get_active_bits(Satellite_mask[0], Satellite_mask_len[0]);
	sat_cnt += get_active_bits(Satellite_mask[1], Satellite_mask_len[1]);
	sig_cnt = get_active_bits(Signal_mask, Signal_mask_len);

	/*	printf("sat_cnt=%d\n", sat_cnt);
	printf("sig_cnt=%d\n", sig_cnt);*/

	/*Cell mask*/
	if (sat_cnt*sig_cnt > 64)
	{
		printf("Invalid CellMask %d\n", sat_cnt*sig_cnt);
		return -1;
	}

	if (k + sat_cnt*sig_cnt > mes_len * 8)
	{
		printf("Invalid message: Len_Cell mask\n");
		return -1;
	}

	offs = 0;
	for (l = j = 0; l < sizeof(Satellite_mask) / sizeof(Satellite_mask[0]); ++l)
	{
		unsigned int mask = Satellite_mask[l];
		unsigned int mask_len = Satellite_mask_len[l];
		for (i = 0; i < mask_len; ++i)
		{
			if (mask & (1 << (mask_len - i - 1)))
			{
				/*printf( "sat[%d]=%d\n", j, i + offs );*/
				sat[j++] = i + offs;
			}
		}
		offs += Satellite_mask_len[l];
	}

	{
		unsigned int mask = Signal_mask;
		unsigned int mask_len = Signal_mask_len;
		for (j = i = 0; i < mask_len; ++i)
		{
			if (mask & (1 << (mask_len - i - 1)))
			{
				/*printf( "sig[%d]=%d\n", j, i );*/
				sig[j++] = i;
			}
		}
	}

	ncell = 0;
	for (i = 0; i < sat_cnt; i++)
	{
		unsigned int mask = getbitu(Raw, k, sig_cnt); k += sig_cnt;
		sat_sig_cnt[sat[i]] = 0;

		for (j = 0; j < sig_cnt; ++j)
		{
			if (mask & (1 << (sig_cnt - j - 1)))
			{
				/*printf( "sat=%d sig[%d]=%d\n", sat[ i], sat_sig_cnt[ sat[ i]], sig[ j] );*/
				sig_type[sat[i]][sat_sig_cnt[sat[i]]] = sig[j];
				sat_sig_cnt[sat[i]]++;
			}
		}

		ncell += sat_sig_cnt[sat[i]];
		/*printf( "sat_sig_cnt[%d]=%d mask=%x\n", sat[ i], sat_sig_cnt[ sat[ i]], mask );*/
	}

	printf("ncell=%d\n", ncell);

	/*Sattelite Data*/
	need_bits = 0;
	if (Nms_follow == 1)
		need_bits += sat_cnt * 8;
	if (Pseudo_range_follow == 2)
		need_bits += sat_cnt * 10;
	if (Supplementary_follow == 2)
		need_bits += sat_cnt * 32;
	if (k + need_bits > mes_len * 8)
	{
		printf("Invalid message: Len_Sattelite Data\n");
		return -1;
	}

	for (i = 0; i < sat_cnt; i++)
	{
		if (Nms_follow == 1)
		{
			Int_num_sat_ranges[sat[i]] = getbitu(Raw, k, 8); k += 8;
		}
		else
			Int_num_sat_ranges[sat[i]] = 255;
	}

	if (Pseudo_range_follow == 2)
	{
		for (i = 0; i < sat_cnt; i++)
		{
			Sat_rough_range[sat[i]] = getbitu(Raw, k, 10); k += 10;
		}
	}

	for (i = 0; i < sat_cnt; i++)  /* Reference PseudoRange */
	{
		Reference_P[sat[i]] = Int_num_sat_ranges[sat[i]] * 1024 + Sat_rough_range[sat[i]];
	}

	if (Supplementary_follow == 2)
	{
		for (i = 0; i < sat_cnt; i++)
		{
			Azimuth[sat[i]] = getbitu(Raw, k, 8); k += 8;
			Elevation[sat[i]] = getbitu(Raw, k, 7); k += 7;
			Doppler[sat[i]] = getbits(Raw, k, 14); k += 14;
			Full_Range_Available[sat[i]] = getbitu(Raw, k, 1); k += 1;
			Satellite_Usage_Status[sat[i]] = getbitu(Raw, k, 2); k += 2;
		}
	}


	/* Signal data*/
	need_bits = 0;
	if (Pseudo_range_follow == 1 || Pseudo_range_follow == 2)
		need_bits += ncell*(Resolution == 0 ? 15 : 20);
	if (Carrier_phase_follow == 2)
		need_bits += ncell*(Resolution == 0 ? 16 : 22);
	if (Carrier_phase_follow == 1 || Carrier_phase_follow == 2)
		need_bits += ncell*(Resolution == 0 ? 8 : 10);
	if (Supplementary_follow == 1 || Supplementary_follow == 2)
		need_bits += ncell*(Resolution == 0 ? 6 : 10);
	if (Supplementary_follow == 2)
		need_bits += ncell*(Resolution == 0 ? 56 : 64);
	if (k + need_bits > mes_len * 8)
	{
		printf("Invalid message: Len_Signal Data\n");
		return -1;
	}

	for (i = 0; i < sat_cnt; i++) /*PseudoRange */
	{
		int s = sat[i];
		int sat_sig = sat_sig_cnt[s];
		for (j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[s][j];
			int len = (Resolution == 0 ? 15 : 20);
			FinePseudoRange[s][ss] = getbitu(Raw, k, len); k += len;
		}
	}

	for (i = 0; i < sat_cnt; i++)
	{
		int s = sat[i];
		int sat_sig = sat_sig_cnt[s];
		for (j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[s][j];
			int len = (Resolution == 0 ? 4 : 10);
			CycSlipCounter[s][ss] = getbitu(Raw, k, len); k += len;
			IntCycPhase[s][ss] = getbitu(Raw, k, 12); k += 12;
		}
	}

	for (i = 0; i < sat_cnt; i++) /*Phase */
	{
		int s = sat[i];
		int sat_sig = sat_sig_cnt[s];
		for (j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[s][j];
			int len = (Resolution == 0 ? 8 : 10);
			FracCycPhase[s][ss] = getbitu(Raw, k, len); k += len;
		}
	}

	for (i = 0; i < sat_cnt; i++)
	{
		int s = sat[i];
		int sat_sig = sat_sig_cnt[s];
		for (j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[s][j];
			int len = (Resolution == 0 ? 6 : 10);
			SNR[s][ss] = getbitu(Raw, k, len); k += len;
		}
	}

	for (i = 0; i < sat_cnt; i++)/* Utochnenii doppler tut */
	{
		int s = sat[i];
		int sat_sig = sat_sig_cnt[s];
		for (j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[s][j];
			int len = (Resolution == 0 ? 56 : 64);
			ChannelNumber[s][ss] = getbitu(Raw, k, 8); k += 8;
			FineDoppler[s][ss] = getbits(Raw, k, 15); k += 15;
			ExtSuppData[s][ss][0] = getbitu(Raw, k, 9); k += 9;
			ExtSuppData[s][ss][1] = getbitu(Raw, k, (len - 32)); k += (len - 32);
		}
	}

	/* Position reference*/
	need_bits = 0;
	if (position_presentation == 1)
		need_bits += 128;
	else
	if (position_presentation == 2)
		need_bits += 152;
	else
	if (position_presentation == 3)
		need_bits += 280;
	if (k + need_bits > mes_len * 8)
	{
		printf("Invalid message: Len_PositionRef Data %d %d %d\n", k, mes_len * 8, k + need_bits);
		return -1;
	}

	k += need_bits;

	printf("\n Message Complete!!! %d %d\n\n", k, mes_len * 8);

	/* rtklib_time */
	int n, prn, s, freq[32], ind[32];
	unsigned char code[32];

	time2str(raw->time, www, 0);
	printf("TIME_before=%s\n", www);
	time2gpst(raw->time, &week); /*calculating week*/
	if (time_tag_extension_type == 0)
	{
		raw->time = gpst2time(week, (day*86400.) + (hour*3600.) + primary_time_tag);
	}
	else
	{
		raw->time = gpst2time(week, primary_time_tag);
	}
	time2str(raw->time, www, 0);
	printf("TIME_after=%s\n", www);

	/* RTKlib fields filling*/
	if (128 == GNSSmask)/*GPS*/
	{
		printf("----------!!!GPS obs message!!!\n");
		for (i = n = 0; i < sat_cnt; i++)
		{
			unsigned int si = sat[i];
			const char* sig_obs[32];
			int index;
			double tt;

			for (m = 0; m < sig_cnt; m++)
			{
				/* signal to rinex obs type */
				sig_obs[m] = msm_sig_gps[sig[m]];
				code[m] = obs2code(sig_obs[m], freq + m);
			}
			sigindex(SYS_GPS, code, freq, sig_cnt, 0, ind);

			prn = si + 1; /* sdvigaem schetchik na 1 */

			if ((s = satno(SYS_GPS, prn))) {
				tt = timediff(raw->obs.data[0].time, raw->time);
				if (fabs(tt) > 1E-9) {
					raw->obs.n = 0;
				}
				index = obsindex(&raw->obs, raw->time, s);
			}


			for (j = 0; j < sat_sig_cnt[si] && j < NFREQ + NEXOBS; j++)
				/*for (j = 0; j < sat_sig_cnt[ si] ; j++) */
			{
				int ss = sig_type[si][j];
				double wl = satwavelen(sat, freq[j] - 1, NULL);
				if (s&&index >= 0 && ind[j] >= 0)
				{
					raw->obs.data[index].LLI[ind[j]] = CycSlipCounter[si][ss];
					raw->obs.data[index].P[ind[j]] = RestorePValue(((double)Reference_P[si])*RANGE_MS / 1024., 655.36, FinePseudoRange[si][ss] * 0.02);
					if (wl > 0.0)
					{
						FinePseudoPhase = (double)(FracCycPhase[si][ss] / (256.) + IntCycPhase[si][ss]);
						raw->obs.data[index].L[ind[j]] = RestorePValue((double)Reference_P[si] * RANGE_MS / (wl*1024.), 4096, FinePseudoPhase);
					}
					if (wl > 0.0)
					{
						raw->obs.data[index].D[ind[j]] = -(Doppler[si] + FineDoppler[si][ss] * 0.0001) / wl;
					}
					raw->obs.data[index].SNR[ind[j]] = SNR[si][ss] / 0.25;
					raw->obs.data[index].code[ind[j]] = code[j];
				}
				n++;
			}
			printf("sat %d added\n", s);
		}
	}
	else if (32 == GNSSmask) /*GLO*/
	{
		printf("----------!!!GLO obs message!!!\n");
		for (i = n = 0; i < sat_cnt; i++)
		{
			unsigned int si = sat[i];
			const char* sig_obs[32];
			int index;
			double tt;

			for (m = 0; m < sig_cnt; m++)
			{
				/* signal to rinex obs type */
				sig_obs[m] = msm_sig_glo[sig[m]];
				code[m] = obs2code(sig_obs[m], freq + m);
			}
			sigindex(SYS_GLO, code, freq, sig_cnt, 0, ind);

			prn = si + 1; /* sdvigaem schetchik na 1 */

			if ((s = satno(SYS_GLO, prn))) {
				tt = timediff(raw->obs.data[0].time, raw->time);
				if (fabs(tt) > 1E-9) {
					raw->obs.n = 0;
				}
				index = obsindex(&raw->obs, raw->time, s);
			}


			for (j = 0; j < sat_sig_cnt[si] && j < NFREQ + NEXOBS; j++)
				/*for (j = 0; j < sat_sig_cnt[ si] ; j++) */
			{
				int ss = sig_type[si][j];
				double wl = satwavelen(sat, freq[j] - 1, NULL);
				if (s&&index >= 0 && ind[j] >= 0)
				{
					raw->obs.data[index].LLI[ind[j]] = CycSlipCounter[si][ss];
					raw->obs.data[index].P[ind[j]] = RestorePValue(((double)Reference_P[si])*RANGE_MS / 1024., 655.36, FinePseudoRange[si][ss] * 0.02);
					if (wl > 0.0)
					{
						FinePseudoPhase = (double)(FracCycPhase[si][ss] / (256.) + IntCycPhase[si][ss]);
						raw->obs.data[index].L[ind[j]] = RestorePValue((double)Reference_P[si] * RANGE_MS / (wl*1024.), 4096, FinePseudoPhase);
					}
					if (wl > 0.0)
					{
						raw->obs.data[index].D[ind[j]] = -(Doppler[si] + FineDoppler[si][ss] * 0.0001) / wl;
					}
					raw->obs.data[index].SNR[ind[j]] = SNR[si][ss] / 0.25;
					raw->obs.data[index].code[ind[j]] = code[j];
				}
				n++;
			}
			printf("sat %d added\n", s);
		}
	}
	/*printf("n=%d, nmax=%d\n", raw->obs.n, raw->obs.nmax);*/

	if (multiple)
		return 0;
	else
		return 1;
}

static int decode_atom_nav(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	/*printf("====== YA V NAV ======\n");*/
	int ret = -1, type = 0;

	unsigned int version, ref_st_ID, NAV_message_type, Standardized_message_number;


	printf("===============================================================NAV\n");

	version = getbitu(Raw, k, 3); k += 3;
	printf("version=%d\n", version);

	ref_st_ID = getbitu(Raw, k, 12); k += 12;
	printf("reference stantion ID=%d\n", ref_st_ID);

	NAV_message_type = getbitu(Raw, k, 9); k += 9;
	printf("NAV message type=%d\n", NAV_message_type);

	if ((NAV_message_type <= 0) || (NAV_message_type >= 3)) /*GPS and GLO eph only*/
		return ret;

	Standardized_message_number = getbitu(Raw, k, 12); k += 12;
	printf("Standardized message number=%d\n", Standardized_message_number);


	switch (Standardized_message_number) {
	case 1019:
		ret = decode_type1019(raw, Raw, k, mes_len); break;
	case 1020:
		ret = decode_type1020(raw, Raw, k, mes_len); break;
	}

	if (ret >= 0) {
		type = Standardized_message_number - 1000;
		if (1 <= type&&type <= 299) raw->msgtype[type]++; else raw->msgtype[0]++;
		printf("!!!!`get nav message`!!!!\n");
		printf("raw->msgtype[%d]=%d\n", type, raw->msgtype[type]);
	}
	return ret;
}

int input_atomf(raw_t *raw, FILE *f)
{
	unsigned char Raw[2048];
	int k, type, ret;
	unsigned int mes_len, mes_num;
	unsigned int mes_sub_num;

start:
	k = 0;
	type = -1; ret = 0;

	while (1)
	{
		if (fread(Raw, 1, 1, f) != 1)
			return -2;
		if (Raw[0] == 0xD3)
			break;
	}

	if (fread(Raw + 1, 1, 2, f) != 2)
		return -2;

	mes_len = getbitu(Raw, 14, 10) + 3;
	if (fread(Raw + 3, 1, mes_len, f) != mes_len)
		return -2;

	if (crc24q(Raw, mes_len) != getbitu(Raw, (mes_len)* 8, 24))
	{
		fseek(f, -mes_len - 2, SEEK_CUR);
		goto start;
	}


	k = 24;
	mes_num = getbitu(Raw, k, 12); k += 12;


	if (mes_num != 4095)
	{
		printf("ne4095\n");
		goto start;
	}

	mes_sub_num = getbitu(Raw, k, 4); k += 4;

	/*decode message*/

	switch (mes_sub_num) {
	case 5:
		ret = decode_atom_nav(raw, Raw, k, mes_len);
		if (-1 >= ret)
			goto start;
		type = ret;
		printf("ret=%d\n", ret);
		break;
	case 7:
		ret = decode_atom_rnx(raw, Raw, k, mes_len);
		if (-1 >= ret)
			goto start;
		type = ret;
		break;
	default: goto start;
	}
	printf("return %d\n", type);
	return type;
}

