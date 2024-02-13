IMPORTANT NOTES AND SUMMARY

This semester we will attempt to recreate the previous semesters findings but with more resolved data

########### Pipeline ##############

1) Jy to beam python (changes y axis)
2) fits files python (trims channels, exclude pix, and makes fits files of isolated peaks)
3) calibrating CH3CN python (uses the rest frequency in the header file to calibrate the velocity)
4) automating casa immoments python (makes 0,1,2,8 moment plots for each region)
5) moment mapping casa fits python (prints moment maps in galactic coords with scale)

############ CASA FUNCTIONS ################

- imregrid (to change from RA Dec to GAL)
- exportfits (changes a casa image file to a fits file)
- imsubimage (can take channels out of a fits file)
- immoments (does moment maps)
