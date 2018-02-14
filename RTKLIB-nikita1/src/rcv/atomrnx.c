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
    int i,j;
    
    for (i=0;i<obs->n;i++) {
        if (obs->data[i].sat==sat) return i; /* field already exists */
    }
    if (i>=MAXOBS) return -1; /* overflow */
    
    /* add new field */
    obs->data[i].time=time;
    obs->data[i].sat=sat;
    for (j=0;j<NFREQ+NEXOBS;j++) {
        obs->data[i].L[j]=obs->data[i].P[j]=0.0;
        obs->data[i].D[j]=0.0;
        obs->data[i].SNR[j]=obs->data[i].LLI[j]=obs->data[i].code[j]=0;
    }
    obs->n++;
    return i;
}

static void sigindex(int sys, const unsigned char *code, const int *freq, int n,
                     const char *opt, int *ind)
{
    int i,nex,pri,pri_h[8]={0},index[8]={0},ex[32]={0},j;
    
    /* test code priority */
	
    for (i=0;i<n;i++) {
        if (!code[i]) continue;
        
        if (freq[i]>NFREQ) { /* save as extended signal if freq > NFREQ */
            ex[i]=1;
            continue;
        }
        /* code priority */
        pri=getcodepri(sys,code[i],opt);
        
        /* select highest priority signal */
        if (pri>pri_h[freq[i]-1]) {
            if (index[freq[i]-1]) ex[index[freq[i]-1]-1]=1;
            pri_h[freq[i]-1]=pri;
            index[freq[i]-1]=i+1;
        }
        else ex[i]=1;
    }
	
    /* signal index in obs data */
    for (i=nex=0;i<n;i++) {
        if (ex[i]==0) ind[i]=freq[i]-1;
        else if (nex<NEXOBS) ind[i]=NFREQ+nex++;
        else { /* no space in obs data */
            trace(2,"rtcm msm: no space in obs data sys=%d code=%d\n",sys,code[i]);
            ind[i]=-1;
        }
#if 0
        trace(2,"sig pos: sys=%d code=%d ex=%d ind=%d\n",sys,code[i],ex[i],ind[i]);
#endif
    }
}

#if 1
static FILE *nik_log;

static int nik_printf(const char *format,...)
{
  char buf[ 256 ];
  void *tbl;
  va_list args;
  unsigned int size;

  va_start( args, format );
  size = vsnprintf( buf, sizeof( buf ), format, args );
  va_end( args );

  if( !nik_log )
	  nik_log = fopen( "nik_log.txt", "wt" );

  if( nik_log )
  {
	  fprintf( nik_log, "%s", buf );
	  fflush( nik_log );
  }

  return 0;
}
#else
static int nik_printf(const char *format,...)
{
  return 0;
}
#endif

#define printf nik_printf

static int get_active_bits( unsigned int mask, int len )
{
	int res;
	for( res = 0; mask > 0 && len > 0 ; mask >>= 1, len-- )
		if(( mask & 1 ))
			res++;
	return( res );
}

#define MAX_SATS 64
#define MAX_SIGS 64
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define RANGE_MS    (CLIGHT*0.001)      /* range in 1 ms */

static int decode_atom_rnx(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	extern const char *msm_sig_gps[32];
    int i, j, l, m;
	unsigned int GNSSmask, version ,ref_st_ID , multiple ,IODS ,smoothing ,position_presentation ,primary_GNSS_system,
	primary_time_tag , time_tag_extension_type,fractional_second ,hour ,day , divergence_free_smoothing_indicator,
	cumulative_session_transmitting_time_indicator ,Data_ID_change_counter,Data_ID_follow,Nms_follow,Supplementary_follow,
	Pseudo_range_follow,Carrier_phase_follow,Resolution ,Reserved1,Reserved2, Satellite_mask[MAX_SATS/32], Signal_mask,
	Satellite_mask_len[MAX_SATS/32],Signal_mask_len,sat_cnt, sig_cnt,sat[MAX_SATS],sat_sig_cnt[MAX_SATS],sig[MAX_SIGS],
	sig_type[MAX_SATS][MAX_SIGS], offs, ncell,need_bits,Int_num_sat_ranges[MAX_SATS],Sat_rough_range[MAX_SATS], Azimuth[MAX_SATS],
	Elevation[MAX_SATS], Doppler[MAX_SATS], Full_Range_Available[MAX_SATS], Satellite_Usage_Status[MAX_SATS], ChannelNumber[MAX_SATS][MAX_SIGS],
	FineDoppler[MAX_SATS][MAX_SIGS], FinePseudoRange[MAX_SATS][MAX_SIGS], CycSlipCounter[MAX_SATS][MAX_SIGS], IntCycPhase[MAX_SATS][MAX_SIGS],
	FracCycPhase[MAX_SATS][MAX_SIGS], SNR[MAX_SATS][MAX_SIGS], ExtSuppData[MAX_SATS][MAX_SIGS][2],Reference_P[MAX_SATS];

	if(k + 58 > mes_len*8)
	{
		printf("Invalid message: Len_Message Header \n");
		return -1;
	}

	printf("===============================================================\n");

	 version = getbitu(Raw, k, 3);k+=3;
	printf("version=%d\n", version); 


	 ref_st_ID = getbitu(Raw, k, 12);k+=12;
	printf("reference stantion ID=%d\n", ref_st_ID); 

	 multiple = getbitu(Raw, k, 1);k+=1;
	printf("multiple message bit=%d\n", multiple); 

	 IODS = getbitu(Raw, k, 3);k+=3;
	printf("IODS=%d\n", IODS); 

	 smoothing = getbitu(Raw, k, 3); k+=3;
	printf("Smoothing interval=%d\n", smoothing);

     position_presentation = getbitu(Raw, k, 2); k+=2;
	printf("Position presentation=%d\n", position_presentation);

	GNSSmask = getbitu(Raw, k, 8);k+=8;
	printf("GNSS mask=%x\n", GNSSmask); 
	if(GNSSmask != 1)
        printf("GNSSmask!=1\n");

	 primary_GNSS_system = getbitu(Raw, k, 3);k+=3;
	printf("Primary GNSS system=%d\n", primary_GNSS_system); 

	 primary_time_tag = getbitu(Raw, k, 12);k+=12;
	printf("Primary time tag=%d\n", primary_time_tag); 

    time_tag_extension_type = getbitu(Raw, k, 1);k+=1;
	printf("Time tag extension type=%d\n", time_tag_extension_type); 

	if(time_tag_extension_type == 1)
	{
		 fractional_second = getbitu(Raw, k, 8);k+=8;
	    printf("Fractional second=%d\n", fractional_second); 
	}
	else 
	{
		 hour = getbitu(Raw, k, 5); k+=5;
		printf("Hour=%d\n", hour); 
		 day = getbitu(Raw, k, 3); k+=3;
		printf("Day=%d\n", day); 
	}

	divergence_free_smoothing_indicator = getbitu(Raw, k, 1); k++;
    printf("Divergence free smoothing indicator=%d\n", divergence_free_smoothing_indicator); 

     cumulative_session_transmitting_time_indicator = getbitu(Raw, k, 7); k+=7;
    printf("Cumulative session transmitting time indicator=%d\n", cumulative_session_transmitting_time_indicator);

    /*Observables mask*/
	if(k + 16 > mes_len*8)
	{
		printf("Invalid message: Len_Observable mask\n");
		return -1;
	}

    Data_ID_change_counter = getbitu(Raw, k, 5); k+=5;
    printf("Data_ID_change_counter=%d\n", Data_ID_change_counter);

    Data_ID_follow = getbitu(Raw, k, 1); k+=1;
	printf("Data_ID_follow=%d\n", Data_ID_follow);

    Nms_follow = getbitu(Raw, k, 1); k+=1;
	printf("Nms_follow=%d\n", Nms_follow);

    Supplementary_follow = getbitu(Raw, k, 2); k+=2;
	printf("Supplementary_follow=%d\n", Supplementary_follow);

    Pseudo_range_follow = getbitu(Raw, k, 2); k+=2;
	printf("Pseudo_range_follow=%d\n", Pseudo_range_follow);

    Carrier_phase_follow = getbitu(Raw, k, 2); k+=2;
	printf("Carrier_phase_follow=%d\n", Carrier_phase_follow);

    Resolution = getbitu(Raw, k, 1); k+=1;
	printf("Resolution=%d\n", Resolution);

    Reserved1 = getbitu(Raw, k, 1); k+=1;
	printf("Reserved1=%d\n", Reserved1);

    Reserved2 = getbitu(Raw, k, 1); k+=1;
	printf("Reserved2=%d\n", Reserved2);

	/*Sat/Sig mask*/
	if(Data_ID_follow != 1)
	{
		printf("Invalid message: Data Id follow\n");
        return -1;
	}

	if(version == 2)
	{
	    if(k + 96 > mes_len*8)
	    {
	    	printf("Invalid message: Len_Sat/Sig MASK1\n");
	    	return -1;
	    }

		Satellite_mask[0] = getbitu(Raw, k, 32); k+=32;
		Satellite_mask[1] = getbitu(Raw, k, 32); k+=32;
		printf("Satellite_mask=%08x%08x\n", Satellite_mask[0], Satellite_mask[1]);
        Satellite_mask_len[0]=32; Satellite_mask_len[1]=32;

        Signal_mask = getbitu(Raw, k, 32); k+=32;
		printf("Signal_mask=%x\n", Signal_mask);
		Signal_mask_len = 32;
	}
	else
	{
	    if(k + 72 > mes_len*8)
	    {
	    	printf("Invalid message: Len_Sat/Sig MASK2\n");
	    	return -1;
	    }

		Satellite_mask[0] = getbitu(Raw, k, 32); k+=32;
		Satellite_mask[1] = getbitu(Raw, k, 8); k+=8;
		printf("Satellite_mask=%08x%02x\n", Satellite_mask[0], Satellite_mask[1]);
		Satellite_mask_len[0]=32; Satellite_mask_len[1]=8;

		Signal_mask = getbitu(Raw, k, 24); k+=24;
		printf("Signal_mask=%x\n", Signal_mask);
		Signal_mask_len = 24;
		k += 8; /* reserved*/
	}

	sat_cnt = get_active_bits( Satellite_mask[0], Satellite_mask_len[0]);
	sat_cnt += get_active_bits( Satellite_mask[1], Satellite_mask_len[1]);
    sig_cnt = get_active_bits( Signal_mask, Signal_mask_len);

	printf("sat_cnt=%d\n", sat_cnt);
	printf("sig_cnt=%d\n", sig_cnt);

	/*Cell mask*/
	if( sat_cnt*sig_cnt > 64)
	{
		printf("Invalid CellMask %d\n", sat_cnt*sig_cnt);
		return -1;
	}

	if(k + sat_cnt*sig_cnt > mes_len*8)
	{
		printf("Invalid message: Len_Cell mask\n");
		return -1;
	}

	offs = 0;
	for( l = j = 0; l < sizeof(Satellite_mask)/sizeof(Satellite_mask[0]); ++l )
	{
		unsigned int mask = Satellite_mask[l];
		unsigned int mask_len = Satellite_mask_len[l];
		for( i = 0; i < mask_len; ++i )
		{
			if( mask & (1<<( mask_len - i - 1 )) )
			{
				printf( "sat[%d]=%d\n", j, i + offs );
				sat[ j++] = i + offs;
			}
		}
		offs += Satellite_mask_len[l];
	}

	{
		unsigned int mask = Signal_mask;
		unsigned int mask_len = Signal_mask_len;
		for( j = i = 0; i < mask_len; ++i )
		{
			if( mask & (1<<( mask_len - i - 1 )) )
			{
				printf( "sig[%d]=%d\n", j, i );
				sig[ j++] = i;
			}
		}
	}

	ncell = 0;
	for(i = 0; i < sat_cnt; i++)
	{
	    unsigned int mask = getbitu(Raw, k, sig_cnt); k+=sig_cnt;
		sat_sig_cnt[ sat[ i]] = 0;

		for( j = 0; j < sig_cnt ; ++j )
		{
			if( mask & (1<<( sig_cnt - j - 1 )))
			{
				printf( "sat=%d sig[%d]=%d\n", sat[ i], sat_sig_cnt[ sat[ i]], sig[ j] );
				sig_type[ sat[ i]][ sat_sig_cnt[ sat[ i]] ] = sig[ j]; 
				sat_sig_cnt[ sat[ i]]++;
			}
		}

		ncell += sat_sig_cnt[ sat[ i]];
		printf( "sat_sig_cnt[%d]=%d mask=%x\n", sat[ i], sat_sig_cnt[ sat[ i]], mask );
	}

	printf( "ncell=%d\n", ncell );

	/*Sattelite Data*/
	need_bits = 0;
	if( Nms_follow == 1 )
	    need_bits += sat_cnt*8;
	if(Pseudo_range_follow == 2)
		need_bits += sat_cnt*10;
	if(Supplementary_follow == 2)
		need_bits += sat_cnt*32;
	if(k + need_bits > mes_len*8)
	{
		printf("Invalid message: Len_Sattelite Data\n");
		return -1;
	}

    for(i = 0; i < sat_cnt; i++)
	{
		if( Nms_follow == 1 )
		{
			Int_num_sat_ranges[sat[i]] = getbitu(Raw, k, 8); k+=8;	
		}
		else
			Int_num_sat_ranges[sat[i]] = 255;
	}

	if(Pseudo_range_follow == 2)
	{
		for(i = 0; i < sat_cnt; i++)
		{
			Sat_rough_range[sat[i]] = getbitu(Raw, k, 10); k+=10;
		}
	}

	for(i = 0; i < sat_cnt; i++)  /* Reference PseudoRange */
	{
		Reference_P[sat[i]] = Int_num_sat_ranges[sat[i]]*1024 + Sat_rough_range[sat[i]];
	}

	if(Supplementary_follow == 2)
	{
		for(i = 0; i < sat_cnt; i++)
		{
		    Azimuth[sat[i]] = getbitu(Raw, k, 8); k+=8;
		    Elevation[sat[i]] = getbitu(Raw, k, 7); k+=7; 
	  	    Doppler[sat[i]] = getbits(Raw, k, 14); k+=14;
	        Full_Range_Available[sat[i]] = getbitu(Raw, k, 1); k+=1;
		    Satellite_Usage_Status[sat[i]] = getbitu(Raw, k, 2); k+=2;
		}
	}


	/* Signal data*/
	need_bits = 0;
	if(Pseudo_range_follow == 1 || Pseudo_range_follow == 2)
		need_bits += ncell*( Resolution == 0 ? 15 : 20);
	if(Carrier_phase_follow == 2)
		need_bits += ncell*( Resolution == 0 ? 16 : 22);
	if(Carrier_phase_follow == 1 || Carrier_phase_follow == 2)
		need_bits += ncell*( Resolution == 0 ? 8 : 10);
	if(Supplementary_follow == 1 || Supplementary_follow == 2)
		need_bits += ncell*( Resolution == 0 ? 6 : 10);
	if(Supplementary_follow == 2)
		need_bits += ncell*( Resolution == 0 ? 56 : 64);
	if(k + need_bits > mes_len*8)
	{
		printf("Invalid message: Len_Signal Data\n");
		return -1;
	}

	for(i = 0; i < sat_cnt; i++) /*PseudoRange */
	{
		int s = sat[ i];
		int sat_sig = sat_sig_cnt[ s];
		for(j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[ s][ j];
			int len = ( Resolution == 0 ? 15 : 20);
			FinePseudoRange[ s][ ss] = getbitu(Raw, k, len); k+=len;
		}
	}

	for(i = 0; i < sat_cnt; i++)
	{
		int s = sat[ i];
		int sat_sig = sat_sig_cnt[ s];
		for(j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[ s][ j];
			int len = ( Resolution == 0 ? 4 : 10);
			CycSlipCounter[ s][ ss] = getbitu(Raw, k, len); k+=len;
			IntCycPhase[ s][ ss] = getbitu(Raw, k, 12); k+=12;
		}
	}

	for(i = 0; i < sat_cnt; i++) /*Phase */
	{
		int s = sat[ i];
		int sat_sig = sat_sig_cnt[ s];
		for(j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[ s][ j];
			int len = ( Resolution == 0 ? 8 : 10);
			FracCycPhase[ s][ ss] = getbitu(Raw, k, len); k+=len;
		}
	}

	for(i = 0; i < sat_cnt; i++)
	{
		int s = sat[ i];
		int sat_sig = sat_sig_cnt[ s];
		for(j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[ s][ j];
			int len = ( Resolution == 0 ? 6 : 10);
			SNR[ s][ ss] = getbitu(Raw, k, len); k+=len;
		}
	}

	for(i = 0; i < sat_cnt; i++)/* Utochnenii doppler tut */
	{
		int s = sat[ i];
		int sat_sig = sat_sig_cnt[ s];
		for(j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[ s][ j];
			int len = ( Resolution == 0 ? 56 : 64);
			ChannelNumber[s][ss] = getbitu(Raw, k, 8); k+=8;
			FineDoppler[s][ss] = getbits(Raw, k, 15); k+=15;
			ExtSuppData[ s][ ss][0] = getbitu(Raw, k, 9); k+=9;
			ExtSuppData[ s][ ss][1] = getbitu(Raw, k, (len-32)); k+=(len-32);
		}
	}

	/* Position reference*/
	need_bits = 0;
	if(position_presentation == 1)
		need_bits += 128;
	else
	if(position_presentation == 2)
		need_bits += 152;
	else
	if(position_presentation == 3)
		need_bits += 280;
	if(k + need_bits > mes_len*8)
	{
		printf("Invalid message: Len_PositionRef Data %d %d %d\n", k, mes_len*8, k + need_bits);
		return -1;
	}

	k += need_bits;

printf("\n Message Complete!!! %d %d\n\n", k, mes_len*8);

    /* rtklib_time */
    int n,prn,s,freq[32],ind[32];
	unsigned char code[32];


	if(time_tag_extension_type == 0)
	{
	    raw->time = gpst2time(0, (day*86400.) + (hour*3600.) + primary_time_tag);
	}
	else
	{
        raw->time = gpst2time(0, primary_time_tag);
	}
    /* RTKlib fields filling*/
    for (i=n=0; i< sat_cnt ;i++) 
	{
		unsigned int si = sat[i];
		const char* sig_obs[32];
		int index;
		double tt;
		
		for(m=0; m < sig_cnt; m++)
		{
			sig_obs[m] = msm_sig_gps[sig[m]];
			code[m]=obs2code(sig_obs[m],freq+m);
		}
		sigindex(SYS_GPS,code,freq,sig_cnt,0,ind);

        prn=si+1; /* sdvigaem schetchik na 1 */

        if ((s=satno(SYS_GPS,prn))) {
            tt=timediff(raw->obs.data[0].time,raw->time);
            if (fabs(tt)>1E-9) {
                raw->obs.n=0;
            }
            index=obsindex(&raw->obs,raw->time,s);
        }

        for (j = 0; j < sat_sig_cnt[ si] && j < NFREQ+NEXOBS; j++) 
		/*for (j = 0; j < sat_sig_cnt[ si] ; j++) */
		{
			int ss = sig_type[ si ][ j];
			double wl=satwavelen(sat,freq[j]-1,NULL);
			if (s&&index>=0&&ind[j]>=0)
			{
				raw->obs.data[index].LLI[ind[j]] = CycSlipCounter[ si][ ss];
				raw->obs.data[index].P[ind[j]] = RestorePValue( ((double)Reference_P[si])*RANGE_MS/1024., 655.36,  FinePseudoRange[ si][ ss]*0.02);
				if(wl > 0.0)
				{
                                    double FinePseudoPhase = FracCycPhase[ si][ ss]/(256.) + IntCycPhase[si][ss];
				    raw->obs.data[index].L[ind[j]] = RestorePValue( ((double)Reference_P[si])*RANGE_MS/(wl*1024.))  , 4096, FinePseudoPhase);
				}
				if(wl > 0.0)
				{
				    raw->obs.data[index].D[ind[j]] = -(Doppler[si] + FineDoppler[si][ss]*0.0001)/wl;
				}
				raw->obs.data[index].SNR[ind[j]] = SNR[ si][ ss]/0.25;
				raw->obs.data[index].code[ind[j]] = code[j];
			}
			n++;
        }

		printf( "sat %d added\n", s );
    }
	return 1;
}

static int decode_atom_nav(raw_t *raw, unsigned char *Raw, int k, int mes_len)
{
	printf("====== YA V NAV ======\n");
	return 1;
}

int input_atomf(raw_t *raw, FILE *f)
{
	unsigned char Raw[ 2048 ];
        int k;
	unsigned int mes_len, mes_num;
	unsigned int mes_sub_num;


start:
    k=0;

    while(1)
	{
        if(fread(  Raw, 1, 1, f ) != 1)
			return -2;
		if(Raw[0] == 0xD3)
          break;
	}

    if(fread(  Raw+1, 1, 2, f ) != 2)
		return -2;

    mes_len = getbitu(Raw, 14, 10) + 3;
	if(fread(  Raw + 3, 1, mes_len , f ) != mes_len)
		return -2;

	if(crc24q(Raw, mes_len) != getbitu(Raw, (mes_len)*8, 24))
	{
		fseek( f, -mes_len-2 , SEEK_CUR );
		goto start;
	}

	k = 24;
	mes_num = getbitu(Raw, k, 12); k+=12;
	if( mes_num != 4095)
	{
		printf("ne4095\n");
		goto start;
	}

	 mes_sub_num = getbitu(Raw, k, 4); k+=4;

	/*decode message*/

	switch(mes_sub_num) {
		case 5: 
			if(-1 == decode_atom_nav(raw, Raw, k, mes_len) )
				goto start;
			break;
		case 7: 
			if(-1 == decode_atom_rnx(raw, Raw, k, mes_len) )
				goto start;
			break;
		default: goto start;
	}

	return 1;
}

