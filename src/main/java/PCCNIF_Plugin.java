import ij.*;
import ij.plugin.*;
import java.io.File;
public class PCCNIF_Plugin implements PlugIn {

    // Paths and file names
    private String workPath = "C:\\Data\\PDI";
    private String resultsPath = workPath + "\\results";
    private String expImageFileName = "Exp061_D5_1_dcxcy3_HAcy2_bIIItubcy5_dapi_20X.czi";

    public void run(String arg) {

        String filename = workPath + "\\" + expImageFileName;

        //Open czi image
        IJ.run("Bio-Formats Importer", "open=" + filename + " autoscale color_mode=Grayscale open_files split_channels view=Hyperstack stack_order=XYCZT");


        // Close channels 1 and 2 since we do not use them
        IJ.selectWindow(expImageFileName + " - C=1");
        IJ.run("Close", "");
        IJ.selectWindow(expImageFileName + " - C=2");
        IJ.run("Close", "");

        //Generate a nuclei mask from the blue channel
        ImagePlus blueChannel = WindowManager.getImage(expImageFileName + " - C=0");
        ImagePlus nucleiMask = generateNucleiMask(blueChannel);

        //apply the nuclei mask to the green channel
        ImagePlus greenChImg = WindowManager.getImage(expImageFileName + " - C=3");
        ImagePlus minImg = applyMaskToGreenChannel(nucleiMask, greenChImg);


        //close all unnecessary windows
        /*
        nucleiMask.changes= false;
        nucleiMask.show();
        IJ.run("Close");

        greenChImg.changes= false;
        greenChImg.show();
        IJ.run("Close");
        */
        //proceed with processing n min image
        minImg.show();

        //Get positive matches
        ImagePlus outlines = countGreenMatches(minImg);


        // Save image calculation
        // create results dir
        File newDir = new File(resultsPath);
        newDir.mkdirs();

        //save image
        IJ.saveAs(outlines, "PNG", resultsPath + "\\" + expImageFileName + " - result.bmp");

        // save count results
        IJ.saveAs("Results", resultsPath + "\\Results.xls");

        //IJ.run("Close", "");
        //outlines.close();


        //TODO measure results characteristics and analyze, classify

    }


    /**
     * Channel is the blue channel that contains the cells' nuclei
     *
     * @return mask a binary mask representing the nuclei morphology for all cells
     */
    private ImagePlus generateNucleiMask(ImagePlus blueChannel) {

        // pre-processing //
        IJ.run(blueChannel, "Enhance Contrast...", "saturated=0.9 equalize");

        // processing //

        //subtract background
        IJ.run(blueChannel, "Subtract Background...", "rolling=50");

        //threshold with black background
        IJ.setAutoThreshold(blueChannel, "Default dark");
        Prefs.blackBackground = true;

        //convert to mask
        IJ.run(blueChannel, "Convert to Mask", "");

        // apply closing operation to fill holes in cell nuclei
        IJ.run(blueChannel, "Close-", "");

        //apply watershed filter to separate nuclei clusters
        IJ.run(blueChannel, "Watershed", "");

        return blueChannel;
    }


    // Apply segmentation mask from blue channel to green channel to eliminate information outside the cells nuclei
    private ImagePlus applyMaskToGreenChannel(ImagePlus mask, ImagePlus greenChImg) {

        IJ.run(greenChImg, "8-bit", "");
        ImageCalculator ic = new ImageCalculator();
        return ic.run("Min create", mask, greenChImg);

    }


    /**
     * @param imp the segmented image of channel 3 with the channel 0 mask
     * @return the outline of all positive matches
     */
    public ImagePlus countGreenMatches(ImagePlus imp) {

        IJ.run(imp, "Enhance Contrast...", "saturated=0.9 equalize");
        IJ.setAutoThreshold(imp, "Default dark");
        Prefs.blackBackground = true;
        IJ.run(imp, "Convert to Mask", "");
        IJ.run(imp, "Analyze Particles...", "size=10-Infinity circularity=0.40-1.00 show=Outlines display clear record in_situ");

        //return imp modified
        return imp;

    }

}