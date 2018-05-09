#include <stdio.h>
#include <stdarg.h>
#include "rtklib.h"

#define P2_10       0.0009765625          /* 2^-10 */
#define P2_34       5.820766091346740E-11 /* 2^-34 */
#define P2_46       1.421085471520200E-14 /* 2^-46 */
#define P2_59       1.734723475976810E-18 /* 2^-59 */

/* adjust weekly rollover of bdt time ----------------------------------------*/
static int adjbdtweek(int week)
{
	int w;
	(void)time2bdt(gpst2bdt(utc2gpst(timeget())), &w);
	if (w<1) w = 1; /* use 2006/1/1 if time is earlier than 2006/1/1 */
	return week + (w - week + 512) / 1024 * 1024;
}

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
		nik_log = fopen("nik_log.nik", "wt");

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

/* get signed 38bit field ----------------------------------------------------*/
static double getbits_38(const unsigned char *buff, int pos)
{
	return (double)getbits(buff, pos, 32)*64.0 + getbitu(buff, pos + 32, 6);
}

/* decode GPS ephemerides -----------------------------------------*/
static int decode_GPS_eph(raw_t *raw, unsigned char *Raw, int k, int mes_len)
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

/* decode GLO ephemerides -------------------------------------*/
static int decode_GLO_eph(raw_t *raw, unsigned char *Raw, int k, int mes_len)
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

/* decode SBAS ephemerides -----------------------------------------*/
static int decode_SBAS_eph(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	return 2;
}

/* decode GAL ephemerides -----------------------------------------*/
static int decode_GAL_eph(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	eph_t eph = { 0 };
	double toc, sqrtA;
	char *msg;
	int i = k, prn, sat, week, e5a_hs, e5a_dvs, rsv, sys = SYS_GAL;

	if (i + 484 <= mes_len * 8) {
		prn = getbitu(Raw, i, 6);              i += 6;
		week = getbitu(Raw, i, 12);              i += 12;
		eph.iode = getbitu(Raw, i, 10);              i += 10;
		eph.sva = getbitu(Raw, i, 8);              i += 8;
		eph.idot = getbits(Raw, i, 14)*P2_43*SC2RAD; i += 14;
		toc = getbitu(Raw, i, 14)*60.0;         i += 14;
		eph.f2 = getbits(Raw, i, 6)*P2_59;        i += 6;
		eph.f1 = getbits(Raw, i, 21)*P2_46;        i += 21;
		eph.f0 = getbits(Raw, i, 31)*P2_34;        i += 31;
		eph.crs = getbits(Raw, i, 16)*P2_5;         i += 16;
		eph.deln = getbits(Raw, i, 16)*P2_43*SC2RAD; i += 16;
		eph.M0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.cuc = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.e = getbitu(Raw, i, 32)*P2_33;        i += 32;
		eph.cus = getbits(Raw, i, 16)*P2_29;        i += 16;
		sqrtA = getbitu(Raw, i, 32)*P2_19;        i += 32;
		eph.toes = getbitu(Raw, i, 14)*60.0;         i += 14;
		eph.cic = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.OMG0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.cis = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.i0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.crc = getbits(Raw, i, 16)*P2_5;         i += 16;
		eph.omg = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.OMGd = getbits(Raw, i, 24)*P2_43*SC2RAD; i += 24;
		eph.tgd[0] = getbits(Raw, i, 10)*P2_32;        i += 10; /* E5a/E1 */
		e5a_hs = getbitu(Raw, i, 2);              i += 2; /* OSHS */
		e5a_dvs = getbitu(Raw, i, 1);              i += 1; /* OSDVS */
		rsv = getbitu(Raw, i, 7);
	}
	else {

		return -1;
	}

	
	if (!(sat = satno(sys, prn))) {
		trace(2, "rtcm3 1045 satellite number error: prn=%d\n", prn);
		return -1;
	}
	eph.sat = sat;
	eph.week = adjgpsweek(week % 1024);
	eph.toe = gpst2time(eph.week, eph.toes);
	eph.toc = gpst2time(eph.week, toc);
	eph.ttr = raw->time;
	eph.A = sqrtA*sqrtA;
	eph.svh = (e5a_hs << 4) + (e5a_dvs << 3);
	eph.code = 2; /* data source = f/nav e5a */

	raw->nav.eph[sat - 1] = eph;
	raw->ephsat = sat;
	return 2;
}

/* decode QZSS ephemerides -----------------------------------------*/
static int decode_QZSS_eph(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	eph_t eph = { 0 };
	double toc, sqrtA;
	char *msg;
	int i = k, prn, sat, week, sys = SYS_QZS;

	if (i + 473 <= mes_len * 8) {
		prn = getbitu(Raw, i, 4) + 192;          i += 4;
		toc = getbitu(Raw, i, 16)*16.0;         i += 16;
		eph.f2 = getbits(Raw, i, 8)*P2_55;        i += 8;
		eph.f1 = getbits(Raw, i, 16)*P2_43;        i += 16;
		eph.f0 = getbits(Raw, i, 22)*P2_31;        i += 22;
		eph.iode = getbitu(Raw, i, 8);              i += 8;
		eph.crs = getbits(Raw, i, 16)*P2_5;         i += 16;
		eph.deln = getbits(Raw, i, 16)*P2_43*SC2RAD; i += 16;
		eph.M0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.cuc = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.e = getbitu(Raw, i, 32)*P2_33;        i += 32;
		eph.cus = getbits(Raw, i, 16)*P2_29;        i += 16;
		sqrtA = getbitu(Raw, i, 32)*P2_19;        i += 32;
		eph.toes = getbitu(Raw, i, 16)*16.0;         i += 16;
		eph.cic = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.OMG0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.cis = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.i0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.crc = getbits(Raw, i, 16)*P2_5;         i += 16;
		eph.omg = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.OMGd = getbits(Raw, i, 24)*P2_43*SC2RAD; i += 24;
		eph.idot = getbits(Raw, i, 14)*P2_43*SC2RAD; i += 14;
		eph.code = getbitu(Raw, i, 2);              i += 2;
		week = getbitu(Raw, i, 10);              i += 10;
		eph.sva = getbitu(Raw, i, 4);              i += 4;
		eph.svh = getbitu(Raw, i, 6);              i += 6;
		eph.tgd[0] = getbits(Raw, i, 8)*P2_31;        i += 8;
		eph.iodc = getbitu(Raw, i, 10);              i += 10;
		eph.fit = getbitu(Raw, i, 1) ? 0.0 : 2.0; /* 0:2hr,1:>2hr */
	}
	else {

		return -1;
	}

	if (!(sat = satno(sys, prn))) {
		trace(2, "rtcm3 1044 satellite number error: prn=%d\n", prn);
		return -1;
	}
	eph.sat = sat;
	eph.week = adjgpsweek(week);
	eph.toe = gpst2time(eph.week, eph.toes);
	eph.toc = gpst2time(eph.week, toc);
	eph.ttr = raw->time;
	eph.A = sqrtA*sqrtA;

		if (eph.iode == raw->nav.eph[sat - 1].iode&&
			eph.iodc == raw->nav.eph[sat - 1].iodc) return 0; /* unchanged */

	raw->nav.eph[sat - 1] = eph;
	raw->ephsat = sat;
	return 2;
}

/* decode BeiDou ephemerides -----------------------------------------*/
static int decode_BeiDou_eph(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	eph_t eph = { 0 };
	double toc, sqrtA;
	char *msg;
	int i = k, prn, sat, week, sys = SYS_CMP;

	if (i + 476 <= mes_len * 8) {
		prn = getbitu(Raw, i, 6);              i += 6;
		week = getbitu(Raw, i, 10);              i += 10;
		eph.sva = getbitu(Raw, i, 4);              i += 4;
		eph.code = getbitu(Raw, i, 2);              i += 2;
		eph.idot = getbits(Raw, i, 14)*P2_43*SC2RAD; i += 14;
		eph.iode = getbitu(Raw, i, 8);              i += 8;
		toc = getbitu(Raw, i, 16)*16.0;         i += 16;
		eph.f2 = getbits(Raw, i, 8)*P2_55;        i += 8;
		eph.f1 = getbits(Raw, i, 16)*P2_43;        i += 16;
		eph.f0 = getbits(Raw, i, 22)*P2_31;        i += 22;
		eph.iodc = getbitu(Raw, i, 10);              i += 10;
		eph.crs = getbits(Raw, i, 16)*P2_5;         i += 16;
		eph.deln = getbits(Raw, i, 16)*P2_43*SC2RAD; i += 16;
		eph.M0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.cuc = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.e = getbitu(Raw, i, 32)*P2_33;        i += 32;
		eph.cus = getbits(Raw, i, 16)*P2_29;        i += 16;
		sqrtA = getbitu(Raw, i, 32)*P2_19;        i += 32;
		eph.toes = getbitu(Raw, i, 16)*16.0;         i += 16;
		eph.cic = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.OMG0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.cis = getbits(Raw, i, 16)*P2_29;        i += 16;
		eph.i0 = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.crc = getbits(Raw, i, 16)*P2_5;         i += 16;
		eph.omg = getbits(Raw, i, 32)*P2_31*SC2RAD; i += 32;
		eph.OMGd = getbits(Raw, i, 24)*P2_43*SC2RAD; i += 24;
		eph.tgd[0] = getbits(Raw, i, 8)*P2_31;        i += 8;
		eph.svh = getbitu(Raw, i, 6);              i += 6;
		eph.flag = getbitu(Raw, i, 1);              i += 1;
		eph.fit = getbitu(Raw, i, 1) ? 0.0 : 4.0; /* 0:4hr,1:>4hr */
	}
	else {

		return -1;
	}

	
	if (!(sat = satno(sys, prn))) {
		trace(2, "rtcm3 1047 satellite number error: prn=%d\n", prn);
		return -1;
	}
	eph.sat = sat;
	eph.week = adjbdtweek(week);
	eph.toe = bdt2gpst(bdt2time(eph.week, eph.toes)); /* bdt -> gpst */
	eph.toc = bdt2gpst(bdt2time(eph.week, toc));      /* bdt -> gpst */
	eph.ttr = raw->time;
	eph.A = sqrtA*sqrtA;

		if (eph.iode == raw->nav.eph[sat - 1].iode) return 0; /* unchanged */
	
	raw->nav.eph[sat - 1] = eph;
	raw->ephsat = sat;
	return 2;
}

/* decode ATOM.RNX message -----------------------------------------*/
static int decode_atom_rnx(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	extern const char *msm_sig_gps[32], *msm_sig_glo[32], *msm_sig_qzs[32], *msm_sig_gal[32], *msm_sig_sbs[32], *msm_sig_cmp[32];
	char www[64];
	int i, j, l, m, week, sys, g=7, gnss[8];
	double FinePseudoPhase, X_coordinate, Y_coordinate, Z_coordinate, Antenna_height,
		X_velocity, Y_velocity, Z_velocity, Receiver_clock_offset, Receiver_clock_drift;
	int Doppler[MAX_SATS], FineDoppler[MAX_SATS][MAX_SIGS];
	unsigned int GNSSmask, version, ref_st_ID, multiple, IODS, smoothing, position_presentation, primary_GNSS_system,
		primary_time_tag, time_tag_extension_type, fractional_second, hour, day, divergence_free_smoothing_indicator,
		cumulative_session_transmitting_time_indicator, Data_ID_change_counter, Data_ID_follow, Nms_follow, Supplementary_follow,
		Pseudo_range_follow, Carrier_phase_follow, Resolution, Reserved1, Reserved2, Satellite_mask[MAX_SATS / 32], Signal_mask,
		Satellite_mask_len[MAX_SATS / 32], Signal_mask_len, sat_cnt, sig_cnt, sat[MAX_SATS], sat_sig_cnt[MAX_SATS], sig[MAX_SIGS],
		sig_type[MAX_SATS][MAX_SIGS], offs, ncell, need_bits, Int_num_sat_ranges[MAX_SATS], Sat_rough_range[MAX_SATS], Azimuth[MAX_SATS],
		Elevation[MAX_SATS], Full_Range_Available[MAX_SATS], Satellite_Usage_Status[MAX_SATS], ChannelNumber[MAX_SATS][MAX_SIGS],
		FinePseudoRange[MAX_SATS][MAX_SIGS], CycSlipCounter[MAX_SATS][MAX_SIGS], IntCycPhase[MAX_SATS][MAX_SIGS],
		FracCycPhase[MAX_SATS][MAX_SIGS], SNR[MAX_SATS][MAX_SIGS], ExtSuppData[MAX_SATS][MAX_SIGS][2], Reference_P[MAX_SATS], 
		Motion_flag, Position_quality_flag, Position_tagging, Clarifier_switch, Clarification_data, Clock_estimate_status, 
		ITRF_epoch_year, GPS_UTC_time_offset, The_number_of_GNSS_time_cycles, Receiver_time_status;

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
	for (j = 1; j < 129; j *= 2)
	{
		gnss[g] = !(!(GNSSmask&j));
		printf("%d", !(!(GNSSmask&j)));
		g--;
	}
	printf("\n");
	if ((0 == GNSSmask) && (0 == position_presentation))
	{
		printf("no info message\n");
		printf("raw->obs.n=%d\n", raw->obs.n);
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
	
	/*Message starts here*/
	for (g = 0; g < 8; g++)
	{
		if (0 == gnss[g])
			continue;
		else
		{
			switch (g){
			case 0: sys = SYS_GPS; break;
			case 1: sys = SYS_SBS; break;
			case 2: sys = SYS_GLO; break;
			case 3: sys = SYS_GAL; break;
			case 4: sys = SYS_QZS; break;
			case 5: sys = SYS_CMP; break;
			default: continue;
			}
		}
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
		printf("sat_cnt=%d\n", sat_cnt);
		printf("sig_cnt=%d\n", sig_cnt);
		
		printf("sat mask=");
		for (j = 0; j < 2; ++j)
		{
			for (i = 0; i < 32; ++i)
			{
				printf("%d", !(!(Satellite_mask[j] & (1 << 32 - i - 1))));
			}
		}
		printf("\n");

		printf("sig mask=");
			for (i = 0; i < 32; ++i)
			{
				printf("%d", !(!(Signal_mask & (1 << 32 - i - 1))));
			}
		printf("\n");


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
					printf( "sat[%d]=%d\n", j, i + offs );
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
		/*printf("CellMask=");*/
		unsigned int *CellMask;
		if (NULL == (CellMask = (unsigned int*)malloc(sat_cnt*sizeof(unsigned int))))
		{
			printf("Error allocating CellMask memory\n");
			return -1;
		}
		for (i = 0; i < sat_cnt; i++)
		{
			CellMask[i] = getbitu(Raw, k, sig_cnt); k += sig_cnt;			
			sat_sig_cnt[sat[i]] = 0;

			for (j = 0; j < sig_cnt; ++j)
			{
				/*printf("%d", !(!(CellMask[i] & (1 << (sig_cnt - j - 1)))) );*/
				if (CellMask[i] & (1 << (sig_cnt - j - 1)))
				{
					/*printf( "sat=%d sig[%d]=%d\n", sat[ i], sat_sig_cnt[ sat[ i]], sig[ j] );*/
					sig_type[sat[i]][sat_sig_cnt[sat[i]]] = sig[j];
					sat_sig_cnt[sat[i]]++;
				}
			}

			ncell += sat_sig_cnt[sat[i]];
			/*printf("sat_sig_cnt[%d]=%d mask=%x\n", sat[i], sat_sig_cnt[sat[i]], CellMask[i]);*/
		}
		printf("\n");


		/*printf("ncell=%d\n", ncell);*/

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

			for (i = n = 0; i < sat_cnt; i++)
			{
				unsigned int si = sat[i];
				const char* sig_obs[32];
				int index;
				double tt, wl;
							
				for (m = 0; m < sig_cnt; m++)
				{
					/* signal to rinex obs type */
					switch (sys) {
					case SYS_GPS: sig_obs[m] = msm_sig_gps[sig[m]]; break;
					case SYS_GLO: sig_obs[m] = msm_sig_glo[sig[m]]; break;
					case SYS_GAL: sig_obs[m] = msm_sig_gal[sig[m]]; break;
					case SYS_QZS: sig_obs[m] = msm_sig_qzs[sig[m]]; break;
					case SYS_SBS: sig_obs[m] = msm_sig_sbs[sig[m]]; break;
					case SYS_CMP: sig_obs[m] = msm_sig_cmp[sig[m]]; break;
					default: sig_obs[m] = ""; break;
					}
					code[m] = obs2code(sig_obs[m], freq + m);
				}

				/* freqency index for beidou */
				if (sys == SYS_CMP) {
					if (freq[m] == 5) freq[m] = 2; /* B2 */
					else if (freq[m] == 4) freq[m] = 3; /* B3 */
				}

				sigindex(sys, code, freq, sig_cnt, 0, ind);

				prn = si + 1; /* sdvigaem schetchik na 1 */
				if (sys == SYS_QZS) prn += MINPRNQZS - 1;
				else if (sys == SYS_SBS) prn += MINPRNSBS - 1;

				if ((s = satno(sys, prn))) {
					tt = timediff(raw->obs.data[0].time, raw->time);
					if (fabs(tt) > 1E-9) {
						raw->obs.n = 0;
					}
					index = obsindex(&raw->obs, raw->time, s);
				}

				printf("sat_sig_cnt[%d]=%d\n", si, sat_sig_cnt[si]);
				/*for (j = 0; j < sat_sig_cnt[si] && j < NFREQ + NEXOBS; j++)*/
				int l = 0;
				for (j = 0; j < sat_sig_cnt[si] && j < NFREQ + NEXOBS; j++)
				{
					if ( (j + l) >= sig_cnt)
						break;
					int ss = sig_type[si][j];
					while( (!(CellMask[i] & (1 << (sig_cnt - j - l - 1))))  )
					{
						++l;
					}
					raw->obs.data[index].code[ind[j]] = code[j+l];

					wl = satwavelen(s, freq[j+l] - 1, &raw->nav);

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
							printf("Doppler[%d]=%lf\n", ind[j], raw->obs.data[index].D[ind[j]]);
						}
						raw->obs.data[index].SNR[ind[j]] = SNR[si][ss] / 0.25;
						printf("%d signal dobavlen. code=%d\n", j, code[j+l]);
					}
					else printf("index=%d, ind[%d]=%d\n", index, j, ind[j]);
					n++;
				}
				printf("sat %d added\n", s);
			}
			free(CellMask);
	}

	/*Position message*/
	if (position_presentation) 
	{
		if (k + 128 > mes_len * 8)
		{
			printf("Invalid message: Len_Reference Position \n");
			return -1;
		}

		Motion_flag = getbitu(Raw, k, 1); k += 1;
		printf("Motion_flag=%d\n", Motion_flag);

		Position_quality_flag = getbitu(Raw, k, 3); k += 3; k += 7;
		printf("Position_quality_flag=%d\n", Position_quality_flag);

		Position_tagging = getbitu(Raw, k, 3); k += 3;
		printf("Position_tagging=%d\n", Position_tagging);

		raw->sta.deltype = 0;

		X_coordinate = getbits_38(Raw, k)*0.0001; k += 38;
		printf("X_coordinate=%lf\n", X_coordinate);
		raw->sta.pos[0] = X_coordinate;
		raw->sta.del[0] = 0.;

		Y_coordinate = getbits_38(Raw, k)*0.0001; k += 38;
		printf("Y_coordinate=%lf\n", Y_coordinate);
		raw->sta.pos[1] = Y_coordinate;
		raw->sta.del[1] = 0.;

		Z_coordinate = getbits_38(Raw, k)*0.0001; k += 38;
		printf("Z_coordinate=%lf\n", Z_coordinate);
		raw->sta.pos[2] = Z_coordinate;
		raw->sta.del[2] = 0.;

		if (1 < position_presentation)
		{
			if (k + 24 > mes_len * 8)
			{
				printf("Invalid message: Len_Reference Position \n");
				return -1;
			}
			Clarifier_switch = getbitu(Raw, k, 2); k += 2;
			printf("Clarifier_switch=%d\n", Clarifier_switch);

			if (0 == Clarifier_switch)
			{
				ITRF_epoch_year = getbitu(Raw, k, 6); k += 6;
				printf("ITRF_epoch_year=%d\n", ITRF_epoch_year);
				raw->sta.itrf = ITRF_epoch_year;

				Antenna_height = getbitu(Raw, k, 16)*0.0001; k += 16;
				printf("Antenna_height=%.16lf\n", Antenna_height);
				raw->sta.hgt = Antenna_height;
			}
			else if (1 == Clarifier_switch)
			{
				GPS_UTC_time_offset = getbitu(Raw, k, 6); k += 6;
				printf("GPS_UTC_time_offset=%d\n", GPS_UTC_time_offset);

				The_number_of_GNSS_time_cycles = getbitu(Raw, k, 12); k += 12;
				printf("The_number_of_GNSS_time_cycles=%d\n", The_number_of_GNSS_time_cycles);

				Receiver_time_status = getbitu(Raw, k, 4); k += 4;
				printf("Receiver_time_status=%d\n", Receiver_time_status);

			}

			if (2 < position_presentation)
			{
				if (k + 128 > mes_len * 8)
				{
					printf("Invalid message: Len_Reference Position \n");
					return -1;
				}
				X_velocity = getbits(Raw, k, 25)*0.0001; k += 25;
				printf("X_velocity=%lf\n", X_velocity);

				Y_velocity = getbits(Raw, k, 25)*0.0001; k += 25;
				printf("Y_velocity=%lf\n", Y_velocity);

				Z_velocity = getbits(Raw, k, 25)*0.0001; k += 25;
				printf("Z_velocity=%lf\n", Z_velocity);

				Clock_estimate_status = getbitu(Raw, k, 1); k += 1;
				printf("Clock_estimate_status=%d\n", Clock_estimate_status);

				Receiver_clock_offset = getbits(Raw, k, 30)*0.001; k += 30;
				printf("Receiver_clock_offset=%lf\n", Receiver_clock_offset);

				Receiver_clock_drift = getbits(Raw, k, 22)*0.001; k += 22;
				printf("Receiver_clock_drift=%lf\n", Receiver_clock_drift);

			}
		}
	}

	if (multiple)
		return 0;
	else
	{
		printf("raw->obs.n=%d\n", raw->obs.n);
		return 1;
	}
}

/* decode ATOM.NAV message -----------------------------------------*/
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

	if ((NAV_message_type == 1) || (NAV_message_type == 2) || (NAV_message_type == 4)) /*GPS,GLO,GAL eph only*/
	{
		Standardized_message_number = getbitu(Raw, k, 12); k += 12;
		printf("Standardized message number=%d\n", Standardized_message_number);
	}


	switch (NAV_message_type) {
	case 1:
		ret = decode_GPS_eph(raw, Raw, k, mes_len); break;
	case 2:
		ret = decode_GLO_eph(raw, Raw, k, mes_len); break;
/*	case 3:
		ret = decode_SBAS_eph(raw, Raw, k, mes_len); break;*/
	case 4:
		ret = decode_GAL_eph(raw, Raw, k, mes_len); break;
/*	case 5:
		ret = decode_QZSS_eph(raw, Raw, k, mes_len); break;
	case 6:
		ret = decode_BeiDou_eph(raw, Raw, k, mes_len); break;*/
	}

	if (ret >= 0) {
		type = NAV_message_type;
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
		return input_atomf(raw, f);
	}

	k = 24;
	mes_num = getbitu(Raw, k, 12); k += 12;

	if (mes_num != 4095)
	{
		printf("ne4095\n");
		return input_atomf(raw, f);
	}

	mes_sub_num = getbitu(Raw, k, 4); k += 4;

	/*decode message*/

	switch (mes_sub_num) {
	case 5:
		ret = decode_atom_nav(raw, Raw, k, mes_len);
		if (-1 >= ret)
		{
			return input_atomf(raw, f);
		}
		type = ret;
		printf("ret=%d\n", ret);
		break;
	case 7:
		ret = decode_atom_rnx(raw, Raw, k, mes_len);
		if (-1 >= ret)
		{
			return input_atomf(raw, f);
		}
		type = ret;
		break;
	default: return input_atomf(raw, f);
	}
	printf("return %d\n", type);
	return type;
}

