
## orphan the git tree

## add plotting examples fig1 and fig2
    move examples to folder

## test if conda and pip install work

## update readme
   # list packages needed
        check minimal set of packages needed
        update requirements.py
    # update links

## better refer to IPCC data used (for bars in fig3 and fig4)

## rename contributor function classes to avoid ugly dictionary mapping

## how to deal with figure 1 and 2?
    we can currently only create them meaningfully if calibrationdata is available

## simplify data save format, maybe csv, so data can be easily read by others.

## add link to primap website.

/////// done

## split get_data in part for calibration and part for projection

## split up data directory

    # input (not for public)
    # calibration (available for public)
    # projection (filled by project() function)

## move plot scripts to functions, so that figures are easily reproduced in ipynb.

## think about way to hint people to download links
    try:
        readcsv(data)
    except FileError:
        print "Missing file, XX download file here."

## create library file with main functions, such as project()
    ## add seed to project.

    ## make project function usable for both single GMT timeseries and ensembles
