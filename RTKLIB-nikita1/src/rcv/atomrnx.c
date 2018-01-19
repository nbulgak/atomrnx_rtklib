#include <stdio.h>
#include <stdarg.h>
#include "rtklib.h"

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

int input_atomrnxf(raw_t *raw, FILE *f)
{
	unsigned char Raw[ 2048 ];
    int i, j, k, l;
	unsigned int mes_len, mes_num , GNSSmask;
	unsigned int mes_sub_num;
	unsigned int version ;
	unsigned int ref_st_ID ;
	unsigned int multiple ;
	unsigned int IODS ;
	unsigned int smoothing ;
    unsigned int position_presentation ;
	unsigned int primary_GNSS_system;
	unsigned int primary_time_tag ;
    unsigned int time_tag_extension_type;
    unsigned int fractional_second ;
	unsigned int hour ;
	unsigned int day ;
    unsigned int divergence_free_smoothing_indicator;
    unsigned int cumulative_session_transmitting_time_indicator ;
	unsigned int Data_ID_change_counter;
    unsigned int Data_ID_follow;
    unsigned int Nms_follow;
    unsigned int Supplementary_follow;
    unsigned int Pseudo_range_follow;
    unsigned int Carrier_phase_follow;
    unsigned int Resolution ;
    unsigned int Reserved1;
    unsigned int Reserved2;
	unsigned int Satellite_mask[MAX_SATS/32];
    unsigned int Signal_mask;
	unsigned int Satellite_mask_len[MAX_SATS/32];
	unsigned int Signal_mask_len;
	unsigned int sat_cnt, sig_cnt;
	unsigned int sat[MAX_SATS];
	unsigned int sat_sig_cnt[MAX_SATS];
	unsigned int sig[MAX_SIGS];
	unsigned int sig_type[MAX_SATS][MAX_SIGS];
	unsigned int offs;
	unsigned int ncell;
	unsigned int need_bits;
	unsigned int Int_num_sat_ranges[MAX_SATS];
	unsigned int Sat_rough_range[MAX_SATS];
	unsigned int Extended_sat_data[MAX_SATS];
	unsigned int FinePseudoRange[MAX_SATS][MAX_SIGS];
	unsigned int CycSlipCounter[MAX_SATS][MAX_SIGS];
	unsigned int IntCycPhase[MAX_SATS][MAX_SIGS];
	unsigned int FracCycPhase[MAX_SATS][MAX_SIGS];
	unsigned int SNR[MAX_SATS][MAX_SIGS];
	unsigned int ExtSuppData[MAX_SATS][MAX_SIGS][2];

start:
	i=0;k=0;

    while(1)
	{
        if(fread(  Raw, 1, 1, f ) != 1)
			return -2;
		if(Raw[0] == 0xD3)
          break;
	}

    if(fread(  Raw+1, 1, 2, f ) != 2)
		return -2;

	/*printf("test=%d\n", getbitu(Raw, 0, 8));*/
    mes_len = getbitu(Raw, 14, 10) + 3;
/*	printf("mes_len=%d\n", mes_len);*/
	if(fread(  Raw + 3, 1, mes_len , f ) != mes_len)
		return -2;

	/*printf("crc=%x\n",crc24q(Raw, mes_len));*/
	/*printf("crc1=%x\n",getbitu(Raw, mes_len*8, 24));*/
	if(crc24q(Raw, mes_len) != getbitu(Raw, (mes_len)*8, 24))
	{
		fseek( f, -mes_len-2 , SEEK_CUR );
		goto start;
	}

/*	printf( "%02x %02x %02x %02x %02x %02x %02x %02x\n",*/ 
/*		Raw[ 0], Raw[ 1], Raw[ 2], Raw[ 3],*/
/*		Raw[ 4], Raw[ 5], Raw[ 6], Raw[ 7]*/
/*	);*/

	k = 24;
	mes_num = getbitu(Raw, k, 12); k+=12;
/*	printf("mes_num=%d\n", mes_num);*/

	if( mes_num != 4095)
	{
		printf("ne4095\n");
		goto start;
	}

	 mes_sub_num=getbitu(Raw, k, 4); k+=4;
/*	printf("mes_sub_num=%d\n", mes_sub_num );*/
	if(mes_sub_num != 7)
		goto start;
	if(k + 58 > mes_len*8)
	{
		printf("Invalid message: Len_Message Header \n");
		goto start;
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
/*	if(GNSSmask != 1)*/
/*        printf("GNSSmask!=1\n");*/

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
		goto start;
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
        goto start;
	}

	if(version == 2)
	{
	    if(k + 96 > mes_len*8)
	    {
	    	printf("Invalid message: Len_Sat/Sig MASK1\n");
	    	goto start;
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
	    	goto start;
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
		goto start;
	}

	if(k + sat_cnt*sig_cnt > mes_len*8)
	{
		printf("Invalid message: Len_Cell mask\n");
		goto start;
	}

	offs = 0;
	for( l = j = 0; l < sizeof(Satellite_mask)/sizeof(Satellite_mask[0]); ++l )
	{
		unsigned int mask = Satellite_mask[l];
		unsigned int mask_len = Satellite_mask_len[l];
		for( i = 0; i < mask_len; ++i )
		{
			if( mask & (1<<( mask_len-i-1 )) )
			{
				printf( "sat[%d]=%d\n", j, i+offs );
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
			if( mask & (1<<( sig_cnt-j-1 )))
			{
				printf( "sat=%d sig[%d]=%d\n", sat[ i], sat_sig_cnt[ sat[ i]], sig[ j] );
				sig_type[ sat[ i]][ sat_sig_cnt[ sat[ i]]++ ] = sig[ j]; 
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
		goto start;
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

	if(Supplementary_follow == 2)
	{
		for(i = 0; i < sat_cnt; i++)
		{
			Extended_sat_data[sat[i]] = getbitu(Raw, k, 32); k+=32;
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
		goto start;
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

	for(i = 0; i < sat_cnt; i++)
	{
		int s = sat[ i];
		int sat_sig = sat_sig_cnt[ s];
		for(j = 0; j < sat_sig; j++)
		{
			int ss = sig_type[ s][ j];
			int len = ( Resolution == 0 ? 56 : 64);
			ExtSuppData[ s][ ss][0] = getbitu(Raw, k, 32); k+=32;
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
		goto start;
	}

	k += need_bits;

printf("\n Message Complete!!! %d %d\n\n", k, mes_len*8);

    {/* to rtklib*/
    int n,prn,s;

	if(time_tag_extension_type == 0)
	{
	    raw->time = gpst2time(0, (day*86400.) + (hour*3600.) + primary_time_tag);
	}
	else
	{
        raw->time = gpst2time(0, primary_time_tag);
	}

    for (i=n=0;i< sat_cnt && n < MAXOBS;i++) {
		unsigned int si = sat[i];

        prn=si+1;
        if (!(s=satno(SYS_GPS,prn))) {
            continue;
        }

        raw->obs.data[n].time=raw->time;
        raw->obs.data[n].sat=s;

		printf( "sat=%d sig=%d lim=%d\n", si, sat_sig_cnt[ si], NFREQ );
        for (j = 0;j< sat_sig_cnt[ si] && j < NFREQ;j++) {
			int ss = sig_type[ si ][ j];
			/*for(i = 0; i < si; i++)*/
			{
				raw->obs.data[n].LLI[j] = CycSlipCounter[ si][ ss];
				raw->obs.data[n].L[j] = IntCycPhase[ si][ ss] + FracCycPhase[ si][ ss]/256.;
				raw->obs.data[n].P[j]= FinePseudoRange[ si][ ss]*0.02;
				raw->obs.data[n].D[j] = 0.0;
				raw->obs.data[n].SNR[j] = SNR[ si][ ss]/0.25;
				raw->obs.data[n].code[j] = CODE_NONE+ss;
			}
        }

		printf( "sat %d added\n", si );
        n++;
    }
    raw->obs.n=n;

    }

	return 1;
}

