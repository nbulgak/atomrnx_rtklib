
add_library(rtklib STATIC
  rtklib.h
  rtkcmn.c
  rtksvr.c
  rtkpos.c
  geoid.c
  solution.c
  lambda.c
  sbas.c
  stream.c
  rcvraw.c
  rtcm.c
  preceph.c
  options.c
  pntpos.c
  ppp.c
  ppp_ar.c
  ephemeris.c
  rinex.c
  ionex.c
  rtcm2.c
  rtcm3.c
  rtcm3e.c
  qzslex.c
  convrnx.c
  rcv/novatel.c
  rcv/ublox.c
  rcv/ss2.c
  rcv/crescent.c
  rcv/skytraq.c
  rcv/gw10.c
  rcv/javad.c
  rcv/nvs.c
  rcv/binex.c
  rcv/rt17.c
  rcv/septentrio.c
  rcv/atom.c
)

target_compile_definitions(rtklib PUBLIC -DNFREQ=3 -DNEXOBS=3 -DENAGLO -DENAGAL -DENAQZS -DENACMP -DENAIRN -DTRACE )
target_include_directories(rtklib PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#win32 specific
target_link_libraries(rtklib Winmm ws2_32)

