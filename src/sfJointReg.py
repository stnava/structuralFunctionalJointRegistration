import ants
rdir = "/Users/stnava/code/structuralFunctionalJointRegistration/"
id = '2001'
# collect image data
# data/LS2001/unprocessed/3T/T1w_MPR1/LS2001_3T_T1w_MPR1_gdc.nii.gz
t1fn = rdir + 'data/LS' + id + "/unprocessed/3T/T1w_MPR1/LS"+id+"_3T_T1w_MPR1_gdc.nii.gz"
# now the bold data
boldfnsL = rdir + "data/LS2001/LS2001fmri/unprocessed/3T/rfMRI_REST1_LR/LS2001_3T_rfMRI_REST1_LR_gdc.nii.gz"
boldfnsR = rdir + "data/LS2001/LS2001fmri/unprocessed/3T/rfMRI_REST1_RL/LS2001_3T_rfMRI_REST1_RL_gdc.nii.gz"
# get the ref data
reffns1 = rdir + 'data/LS2001/LS2001fmri/unprocessed/3T/rfMRI_REST1_LR/LS2001_3T_rfMRI_REST1_LR_SBRef_gdc.nii.gz'
reffns2 = rdir + 'data/LS2001/LS2001fmri/unprocessed/3T/rfMRI_REST1_RL/LS2001_3T_rfMRI_REST1_RL_SBRef_gdc.nii.gz'

##  Distortion correction (without a field map)

i1 = ants.image_read( reffns1 )
i2 = ants.image_read( reffns2 )
und = ants.build_template( i1, ( i1, i2 ), 3 )
t1 = ants.image_read( t1fn ).n3_bias_field_correction( 8 ).n3_bias_field_correction( 4 )
bmask = ants.get_mask( und, low_thresh = und.mean() * 0.75, high_thresh=1e9, cleanup = 3 ).iMath_fill_holes()
ants.plot( und, bmask, axis=2, overlay_alpha = 0.33 )
# this is a fragile approach - not really recommended - but it is quick
t1mask = ants.get_mask( t1, low_thresh = t1.mean() * 1.1, high_thresh=1e9, cleanup = 5 ).iMath_fill_holes()
ants.plot( t1, t1mask, axis=2, overlay_alpha = 0.33 )
t1rig = ants.registration( und * bmask, t1 * t1mask, "BOLDRigid" )
t1reg = ants.registration( und * bmask, t1 * t1mask, "SyNOnly",
  initialTransform = t1rig['fwdtransforms'],
  synMetric = 'CC', synSampling = 2, regIterations = (5) )
###########################


# The distortion to the T1 is greatly reduced.

# Brain masking
# Use the BOLD mask to extract the brain from the t1

t1maskFromBold = ants.apply_transforms( t1, bmask, t1reg['invtransforms'],
                                      interpolator = 'nearestNeighbor' )
t1 = ants.n4_bias_field_correction( t1, t1mask, 8 ).n4_bias_field_correction( t1mask, 4 )
bmask =  ants.apply_transforms( und, t1mask, t1reg['fwdtransforms'],
  interpolator = 'nearestNeighbor' ).morphology("close",3)
ofn = rdir + "features/LS" + id + "_mask_py.nii.gz"
ants.image_write( bmask, ofn )
t1toBOLD = ants.apply_transforms( und, t1, t1reg['fwdtransforms'] )
ofn = rdir + "features/LS" + id + "_t1ToBold_py.nii.gz"
ants.image_write( t1toBOLD, ofn )

## Tissue segmentation
# a simple method
################################################
qt1 = ants.iMath_truncate_intensity( t1, 0, 0.95, t1mask )
  t1seg = kmeansSegmentation( qt1, 3, t1mask, 0.2 )
  volumes = labelStats( t1seg$segmentation, t1seg$segmentation )
  rownames( volumes ) = c("background",'csf', 'gm', 'wm' )
  volumes$NormVolume = volumes$Volume / sum( volumes$Volume[-1])
  pander::pander( volumes[ , c("LabelValue","Volume","NormVolume")] )

  # if we look and realize this is not good - fix & redo - caused by local hyperintensity
  if ( volumes["wm","NormVolume"] < 0.2 ) {
    t1mask2 = t1mask * thresholdImage( t1seg$segmentation, 1, 2 )
    t1seg = kmeansSegmentation( t1, 3, t1mask2, 0.1 )
    }
}
plot( t1, t1seg$segmentation, axis = 3, window.overlay=c(0,3) )
boldseg = antsApplyTransforms( und, t1seg$segmentation, t1reg$fwdtransforms,
                               interpolator = 'nearestNeighbor' )
plot(und,boldseg,axis=3,alpha=0.5)
########################################
```


## Template mapping

* include prior information e.g. from meta-analysis or anatomy

```{r totemplate}
if ( ! exists( "treg" ) ) {
  data( "powers_areal_mni_itk" )
  myvoxes = 1:nrow( powers_areal_mni_itk )
  anat = powers_areal_mni_itk$Anatomy[myvoxes]
  syst = powers_areal_mni_itk$SystemName[myvoxes]
  Brod = powers_areal_mni_itk$Brodmann[myvoxes]
  xAAL  = powers_areal_mni_itk$AAL[myvoxes]
  if ( ! exists( "ch2" ) )
    ch2 = ants.image_read( getANTsRData( "ch2" ) )
  treg = ants.registration( t1 * t1mask, ch2, 'SyN' )
}
concatx2 = c( treg$invtransforms, t1reg$invtransforms )
pts2bold = antsApplyTransformsToPoints( 3, powers_areal_mni_itk, concatx2,
                                        whichtoinvert = c( T, F, T, F ) )
ptImg = makePointsImage( pts2bold, bmask, radius = 3 )
plot( und, ptImg, axis=3, window.overlay = range( ptImg ) )
bold2ch2 = antsApplyTransforms( ch2, und,  concatx2,
      whichtoinvert = c( T, F, T, F ) )
# check the composite output
plot( bold2ch2 , axis=3 )
###########################################################
```


# Extracting canonical functional network maps


## preprocessing

back to the fmri ...

    * undistort
    * motion correction

then we will be ready for first-level analysis ...

```{r motion}
if ( ! exists( "motcorr" ) ) {
  bold = ants.image_read( boldfnsR )
  avgBold = getAverageOfTimeSeries( bold )
  # map to und
  boldUndTX = ants.registration( und, avgBold, "SyN", regIterations = c(20,10),
                                 synMetric = "CC", synSampling = 2, verbose = F )
  boldUndTS = antsApplyTransforms( und, bold, boldUndTX$fwd, imagetype = 3  )
  avgBold = getAverageOfTimeSeries( boldUndTS )
  motcorr = antsrMotionCalculation( boldUndTS, verbose = F, typeofTransform = 'Rigid' )

  boldL = ants.image_read( boldfnsL )
  avgBoldL = getAverageOfTimeSeries( boldL )
  boldUndTXL = ants.registration( und, avgBoldL, "SyN", regIterations = c(20,10),
                                 synMetric = "CC", synSampling = 2, verbose = F )
  boldUndTSL = antsApplyTransforms( und, boldL, boldUndTXL$fwd, imagetype = 3, verbose=TRUE  )
  avgBoldL = getAverageOfTimeSeries( boldUndTSL )
  motcorrL = antsrMotionCalculation( boldUndTSL, fixed=avgBold, mask = motcorr$moco_mask,
    verbose = T, typeofTransform = 'Rigid' )

  # bind the two runs
  # newmo = iBind( motcorr$moco_img, motcorrL$moco_img, along = 4 )
  }
plot( ts( motcorr$fd$MeanDisplacement ) )
print( names( motcorr ) )
````

### trim the bold time series for "good" time points

* trim the first $k$ time points

* trim the high motion time points and their $k$ neighbors

We can then perform regression only in the "good" time points.

Or we can impute / interpolate ( see the RestingBold article/vignette in ANTsR documentation ).

## network modeling

now we can do some additional level-one analysis to extract relevant networks.

    * default mode
    * motor

```{r denoise}
# use tissue segmentation to guide compcor
# FIXME - provide reference for this approach

getNetworkBeta <-function( motcorrIn, networkName = 'Default Mode' ) {

  csfAndWM = thresholdImage( boldseg, 1, 1 ) +
           ( thresholdImage( boldseg, 3, 3 ) %>% iMath("ME",1) )
  ccmat = timeseries2matrix( motcorrIn$moco_img, csfAndWM )
  noiseu = compcor( ccmat, ncompcor = 10 )
  smth = c( 1.0, 1.0, 1.0, 2.0 ) # this is for sigmaInPhysicalCoordinates = F
  simg = smoothImage( motcorrIn$moco_img, smth, sigmaInPhysicalCoordinates = F )
  gmmat = timeseries2matrix( simg, thresholdImage( boldseg, 2, 2 ) )
  tr = antsGetSpacing( bold )[4]
  gmmatf = frequencyFilterfMRI( gmmat, tr = tr, freqLo = 0.01, freqHi = 0.1,  opt = "trig" )

  goodtimes = rep( TRUE, dim( motcorrIn$moco_img )[4] )
  goodtimes[ 1:10 ] = FALSE  # signal stabilization
  highMotionTimes = which( motcorrIn$fd$MeanDisplacement >= 0.5 )
  for ( highs in -2:2 )
    highMotionTimes = c( highMotionTimes, highMotionTimes + highs )
  highMotionTimes = sort( unique( highMotionTimes ))
  goodtimes[ highMotionTimes ] = FALSE

  dfnpts = which( powers_areal_mni_itk$SystemName == networkName )
  print( paste(networkName, length(dfnpts)))
  dfnmask = maskImage( ptImg, ptImg, level = dfnpts, binarize = T )
  dfnmat = timeseries2matrix( simg, dfnmask )
  dfnmatf = frequencyFilterfMRI( dfnmat, tr = tr, freqLo = 0.01, freqHi = 0.1,  opt = "trig" )
  dfnsignal = rowMeans( data.matrix( dfnmatf ) )
  locmat = data.matrix( data.frame( gmmatf ))[goodtimes,]
  dnz = aslDenoiseR( locmat, dfnsignal[goodtimes], covariates = motcorrIn$fd[goodtimes,],
                     crossvalidationgroups = 8, selectionthresh = 0.2  )
  dfndf = data.frame( signal = dfnsignal, noiseu = noiseu, fd = motcorrIn$fd )
  dfndf = data.frame( signal = dfnsignal[goodtimes], noiseu = dnz$noiseu, fd = dnz$covariates )
  myform = paste( " mat ~ ", paste( names( dfndf ), collapse = "+") )
  dfnmdl = ilr( dfndf, list( mat = locmat ),  myform )
  dfnBetaImg = makeImage(  thresholdImage( boldseg, 2, 2 ),  dfnmdl$tValue["signal",] )
  return( dfnBetaImg )
  }
dfnBetaImgR = getNetworkBeta( motcorr, 'Default Mode' )
dfnBetaImgL = getNetworkBeta( motcorrL, 'Default Mode' )
handR = getNetworkBeta( motcorr, "Sensory/Somatomotor Hand" )
handL = getNetworkBeta( motcorrL, "Sensory/Somatomotor Hand" )
mouthR = getNetworkBeta( motcorr, "Sensory/Somatomotor Mouth" )
mouthL = getNetworkBeta( motcorrL, "Sensory/Somatomotor Mouth" )
mouth = mouthL + mouthR
hand  = handL + handR
# vizBetaImg = getNetworkBeta( "Visual" )
# salBetaImg = getNetworkBeta( "Salience" )
####################################
##############
```

now we can look at the full continuous beta map for the default mode network

```{r dfnBetas}
plot( und, hand, alpha = 1.0, axis = 3, window.overlay = c(10, max(hand )),
      nslices = 24, ncolumns = 6 )
plot( und, mouth, alpha = 1.0, axis = 3, window.overlay = c(10, max(salBetaImg )),
      nslices = 24, ncolumns = 6 )
plot( und, dfnBetaImgR, alpha = 1.0, axis = 3, window.overlay = c(10, max(dfnBetaImg )),
      nslices = 24, ncolumns = 6 )
ofn = paste0( rdir, "features/LS", id, '_vizBetaImg.nii.gz' )
antsImageWrite( vizBetaImg, ofn )
ofn = paste0( rdir, "features/LS", id, '_dfnBetaImg.nii.gz' )
antsImageWrite( dfnBetaImg, ofn )
ofn = paste0( rdir, "features/LS", id, '_undistort.nii.gz' )
antsImageWrite( und, ofn )
ofn = paste0( rdir, "features/LS", id, '_t1ToBold.nii.gz' )
antsImageWrite( t1toBOLD, ofn )
ofn = paste0( rdir, "features/LS", id, '_mask.nii.gz' )
antsImageWrite( bmask, ofn )
```


Now repeat all of this for the next subject by changing the ID.


# Structural functional joint registration

Finally, use all of this data to do a joint registration using both structural
and functional features.

```{r poorMansHyperalignment}
# ants.registration with multivariateExtras
id1 = '2001'
s1f1 = ants.image_read( paste0( rdir, "features/LS", id1, '_dfnBetaImg.nii.gz' ) )
s1f2 = ants.image_read( paste0( rdir, "features/LS", id1, '_undistort.nii.gz' ) )
s1f3 = ants.image_read( paste0( rdir, "features/LS", id1, '_t1ToBold.nii.gz' ) )
s1fv = ants.image_read( paste0( rdir, "features/LS", id1, '_vizBetaImg.nii.gz' ) )
s1mask = ants.image_read( paste0( rdir, "features/LS", id1, '_mask.nii.gz' ) )
id2 = '3026'
s2f1 = ants.image_read( paste0( rdir, "features/LS", id2, '_dfnBetaImg.nii.gz' ) )
s2f2 = ants.image_read( paste0( rdir, "features/LS", id2, '_undistort.nii.gz' ) )
s2f3 = ants.image_read( paste0( rdir, "features/LS", id2, '_t1ToBold.nii.gz' ) )
s2fv = ants.image_read( paste0( rdir, "features/LS", id2, '_vizBetaImg.nii.gz' ) )
s2mask = ants.image_read( paste0( rdir, "features/LS", id2, '_mask.nii.gz' ) )
#############
jrig = ants.registration( s1f3 * s1mask, s2f3  * s2mask, "Affine"  )
ureg = ants.registration( s1f3, s2f3, "SyNOnly", initialTransform = jrig$fwd,  mask = s1mask )
jreg = ants.registration( s1f3, s2f3, "SyNOnly", initialTransform = jrig$fwd,  mask = s1mask,
  multivariateExtras = list(
   list( "mattes", s1fv, s2fv, 1, 32 ),
   list( "mattes", s1f1, s2f1, 1, 32 ) ), verbose = FALSE )
#############
vizWarped = antsApplyTransforms( s1f1, s2fv, ureg$fwdtransforms )
metric = antsrMetricCreate( s1fv, vizWarped, type="Correlation" )
strWarped = antsApplyTransforms( s1f1, s2f3, ureg$fwdtransforms )
smetric = antsrMetricCreate( s1f3, strWarped, type="Correlation" )
print( paste("univar-dfn", antsrMetricGetValue( metric ), 'str', antsrMetricGetValue( smetric ) ) )
#############
dfnWarped = antsApplyTransforms( s1f1, s2f1, jreg$fwdtransforms )
vizWarped = antsApplyTransforms( s1f1, s2fv, jreg$fwdtransforms )
metric = antsrMetricCreate( s1fv, vizWarped, type="Correlation" )
strWarped = antsApplyTransforms( s1f1, s2f3, jreg$fwdtransforms )
smetric = antsrMetricCreate( s1f3, strWarped, type="Correlation" )
print( paste("mulvar-dfn", antsrMetricGetValue( metric ), 'str', antsrMetricGetValue( smetric ) ) )
#############
```

Result:  *the joint registration performs better on both structural and functional metrics*
as measured by correlation of the features (not ideal).

Visualize the "fixed" subject and DFN first.

Then show the "deformed" subject and DFN second.

```{r regviz}
plot( s1f3, s1f1, alpha = 1.0, axis = 3, window.overlay = c(3, max(dfnBetaImg )),
      nslices = 24, ncolumns = 6, doCropping=FALSE )
plot( strWarped, dfnWarped, alpha = 1.0, axis = 3, window.overlay = c(3, max(dfnBetaImg )),
      nslices = 24, ncolumns = 6, doCropping=FALSE  )
```
