printf("Loading cw_analysis.il\n")
; cwevm


/*****************************************************************
 * TEST FUNCTIONS
 *****************************************************************/

defun( cwevm2_test ()
  pavs = vector(-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1)
  vref = vector(-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0,0,0)
  vmeas_mag   = vector(0.388208,0.435562,0.488691,0.548302,0.615199,0.690284,0.774567,0.86921,0.975524,1.09501,1.22938,1.38057,1.55077,1.74244,1.95833,2.20151,2.4754,2.78388,3.13147,3.52383,3.96812,4.47208,5.04201,5.68003,6.36363,7.01236,7.57652,8.03131,8.42095,8.6897,8.83684,8.93731)
  vmeas_phase = vector(75.5874,75.566,75.547,75.5265,75.494,75.4446,75.3933,75.3271,75.2512,75.1623,75.0596,74.9533,74.8485,74.7478,74.6585,74.5885,74.5431,74.5162,74.5071,74.5075,74.4968,74.4426,74.3026,74.0644,73.7472,73.2589,72.5501,71.3563,69.5733,68.0518,67.122,66.336)

  cwevm_core(pavs, vref, vmeas_mag, vmeas_phase)
)

defun( polar (mag angle)
  anglerad = angle*3.14159/180
  vec = complex( mag*cos(anglerad) mag*sin(anglerad))
)

/*****************************************************************
 * SUPPORT FUNCTIONS
 *****************************************************************/

defun( fixSinglePoint (wave)
  ; Check if family first
  if( famIsFamily( wave ) then
    ; Next check to see if unitary length
    ;if( (length(famGetSweepValues( wave )) == 1)
      wave = if(axlIsICRPProcess() wave famValue(wave car(famGetSweepValues( wave ))) )
    ;)
  )
  wave
)

defun( cwevm_1d (wave)
  printf("cwevm()  ------------------------------------\n")
  ;wave = fixSinglePoint(wave)

  xdata = drGetWaveformXVec(wave)
  ydata = drGetWaveformYVec(wave)

  xlen = drVectorLength(xdata)
  ylen = drVectorLength(ydata)

  printf("xlen = %L\tylen = %L\n", xlen, ylen)

  if( (xlen != ylen) then
    printf("Vectors not the same length!\n")
  )

  pavs        = makeVector(xlen 0)
  vmeas_mag   = makeVector(xlen 0)
  vmeas_phase = makeVector(xlen 0)

  for( k1 0 (xlen-1)
    pavs[k1]        = drGetElem(xdata k1)
    vmeas_mag[k1]   = mag(drGetElem(ydata k1))
    vmeas_phase[k1] = phaseDegUnwrapped(drGetElem(ydata k1))
    ;printf("pavs: %L\t\tmag: %L\t phase: %L\n", pavs[k1], vmeas_mag[k1], vmeas_phase[k1])
  )

  evm = cwevm_core(pavs, pavs, vmeas_mag, vmeas_phase)

  ydata     = drCreateVec('double xlen)
  for( k1 0 (xlen-1)
    drSetElem(ydata k1 evm[k1])
    ;drSetElem(ydata k1 drGetElem(xdata k1))
  )

  newwave = drCreateEmptyWaveform()
  drPutWaveformXVec(newwave xdata)
  drPutWaveformYVec(newwave ydata)
  newwave
) ; end cwevm_1d

defun( cwevm (wave)
  if( famIsFamily(wave) then
    result = famMap( 'cwevm_1d wave )
  else
    result = cwevm_1d(wave)
  )
  result
)
  

defun( cwevm_core (pavs, vref, vmeas_mag, vmeas_phase)
  printf("cwevm_core()\n")
  debug = nil
  printf("exec: cwevm2()\n")
  if( debug then
    printf("Dim 1: %d\n", length(pavs))
  )

  pdf_start = -20
  pdf_stop  = 9
  pdf_step  = 1
  max_clip  = 3
  steps_per_db = 1/pdf_step;

  pdf_length      = (pdf_stop-pdf_start+1)/steps_per_db
  pavs_length     = length(pavs)
  num_pavs_valid = length(pavs) - (pdf_stop/pdf_step) + (max_clip/pdf_step) + floor(steps_per_db)
 
  if( debug then
    printf("pavs_length: %L\n", pavs_length)
  )

  rdb = makeVector(pdf_length 0)
  pdf = makeVector(pdf_length 0)
  for(k1 0 (pdf_length-1)
    rdb[k1] = pdf_start + k1*pdf_step
    rlin    = pow(10, rdb[k1]/10.0)
    drlin   = rlin * pow(10, pdf_step/10.0) - rlin;
    pdf[k1] = drlin * rlin * exp(-rlin)
    if( debug then
      printf("k1: %L\trdb: %L\trlin: %L\tdrlin: %L\tpdf: %L\n", k1, rdb[k1], rlin, drlin, pdf[k1])
    )
  )

  ; Get Pin vector in list form so that min and max can be had. 
  lPavs=vectorToList(pavs)
  minPavs=apply('min lPavs)
  maxPavs=apply('max lPavs)
 
  if( debug then
    printf("[1d] min(Pavs): %L\n", minPavs)
    printf("[1d] max(Pavs): %L\n", maxPavs)
  )

  ; Calculate complex gain
  gain_mag   = makeVector(pavs_length 0)
  gain_phase = makeVector(pavs_length 0)
  if( debug then
    printf("len_pavs: %L\tlen_pdf: %L\tlen_gain: %L\n", pavs_length, pdf_length, length(gain_mag))
  )
  for(k1 0 (pavs_length-1)
    vref_ideal     = sqrt(pow(10, (pavs[k1]-30)/10.0)*2*50)
    gain_mag[k1]   = 20*log10( vmeas_mag[k1] / vref_ideal)
    gain_phase[k1] = vmeas_phase[k1]
  )

  evm      = makeVector(pavs_length 1)  ; Initialize with 0 dB evm
  evm_amam = makeVector(pavs_length 0)
  evm_ampm = makeVector(pavs_length 0)

  for( k1 0 (num_pavs_valid-1)
    evm_amam[k1] = 0
    evm_ampm[k1] = 0
    for( k2 0 (pdf_length-1)
      rpavs = pavs[k1] + rdb[k2]
      if( (rpavs >= minPavs ) && ( rpavs <= maxPavs )
        then
          k1_offset = k1 + rdb[k2]*steps_per_db

          ;printf("k1: %L\tk2: %L\trdb: %L\tsteps_per_db: %L\tk1_offset: %L\n", k1, k2, rdb[k2], steps_per_db, k1_offset)
          ;printf("gain_mag[0]: %L\n", gain_mag[0])
          ;printf("gain_mag[k1_offset]: %L\n", gain_mag[k1_offset])
          ;printf("gain_mag[k1]: %L\n", gain_mag[k1])

          delta_mag   = gain_mag[k1_offset]   - gain_mag[k1]
          delta_phase = gain_phase[k1_offset] - gain_phase[k1]
          evm_amam[k1] = evm_amam[k1] + pdf[k2] * pow(delta_mag,2)
          evm_ampm[k1] = evm_ampm[k1] + pdf[k2] * pow(delta_phase,2)
          ;printf("k1: %L\tk2: %L\tpdf: %L\tmag: %L\tang: %L\tevmam: %L\tevmpm: %L\n", k1, k2, pdf[k2], delta_mag, delta_phase, evm_amam[k1], evm_ampm[k1])
      )
    )
    evm_amam[k1] = ln(10)/20.0      * sqrt(evm_amam[k1])
    evm_ampm[k1] = 3.14159/180.0    * sqrt(evm_ampm[k1])
    evm[k1]      = sqrt(pow(evm_amam[k1],2) + pow(evm_ampm[k1],2))

    ; More efficient
    ;evm_amam[k1] = pow(ln(10)/20.0,2)      * evm_amam[k1]
    ;evm_ampm[k1] = pow(3.14159/180.0,2)    * evm_ampm[k1]
    ;evm[k1]      = sqrt(evm_amam[k1] + evm_ampm[k1])

    if( debug then
      printf("pavs: %L\tevmAM: %L\tevmPM: %L\tevmRMS: %L\n",pavs[k1], evm_amam[k1], evm_ampm[k1], evm[k1]) 
    )
  )
  evm
) ; cwevm2

/*****************************************************************
 * DEBUG FUNCTIONS
 *****************************************************************/


defun( neato2 (wave)
  printf("NEATO2:  ------------------------------------\n")
  sweepName = famGetSweepName( wave )
  printf(strcat("Sweep Name:" sweepName "\n") )
  sweepVals = famGetSweepValues( wave )
  printf("Family Length: %L\n" length(sweepVals))

  xdata = drGetWaveformXVec(wave)
  ydata = drGetWaveformYVec(wave)

  xlen = drVectorLength(xdata)
  ylen = drVectorLength(ydata)

  printf("xlen = %L\tylen = %L\n", xlen, ylen)

  wave
)

defun( returnnumber ()
  1.23
)