/*
 * #%L
 * Bio-Formats Plugins for ImageJ: a collection of ImageJ plugins including the
 * Bio-Formats Importer, Bio-Formats Exporter, Bio-Formats Macro Extensions,
 * Data Browser and Stack Slicer.
 * %%
 * Copyright (C) 2006 - 2015 Open Microscopy Environment:
 *   - Board of Regents of the University of Wisconsin-Madison
 *   - Glencoe Software, Inc.
 *   - University of Dundee
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */

import ij.*;
import ij.io.OpenDialog;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import loci.common.Region;
import loci.formats.FormatException;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import java.io.File;
import java.io.IOException;

/**
 * A very simple example of using Bio-Formats in an ImageJ plugin.
 */
public class PCCNIF_Plugin implements PlugIn {

    // Paths and file names
    private String dir = "C:\\Data\\PDI\\";
    private String resultsPath = dir + "results";
    private String name = "Exp035_D5_a_dcxcy3_flagcy2_bIIItubcy5_dapi_20X.czi";
    String id = dir + name;

    public void run(String arg) {


        // Open czi image
        try {
            ImporterOptions options = new ImporterOptions();
            options.setId(id);
            options.setAutoscale(true);
            options.setStitchTiles(true);
            options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT);
            options.setSplitChannels(true);
            options.isViewHyperstack();
            options.setStackOrder("XYCZT");
            ImagePlus[] imps = BF.openImagePlus(options);
            for (ImagePlus imp : imps) imp.show();
        } catch (FormatException exc) {
            IJ.error("Sorry, an error occurred: " + exc.getMessage());
        } catch (IOException exc) {
            IJ.error("Sorry, an error occurred: " + exc.getMessage());
        }

        // Close channels 2 and 3 since we do not use them
        IJ.selectWindow(name + " - C=2");
        IJ.run("Close");
        IJ.selectWindow(name + " - C=3");
        IJ.run("Close");

        // Generate a nuclei mask from the red channel
        ImagePlus redChannel = WindowManager.getImage(name + " - C=0");
        ImagePlus nucleiMask = generateNucleiMask(redChannel);

        // Apply the nuclei mask to the green channel
        ImagePlus greenChImg = WindowManager.getImage(name + " - C=1");
        ImagePlus minImg = applyMaskToGreenChannel(nucleiMask, greenChImg);

        // Close all unnecessary windows
        nucleiMask.changes= false;
        nucleiMask.show();
        IJ.run("Close");

        greenChImg.changes= false;
        greenChImg.show();
        IJ.run("Close");

        // Proceed with processing n min image
        minImg.show();

        // Get positive matches
        ImagePlus outlines = countGreenMatches(minImg);


        // Save image calculation //

        // Create results dir
        File newDir = new File(resultsPath);
        newDir.mkdirs();

        // Save image
        IJ.saveAs(outlines, "PNG", resultsPath + "\\" + name + " - result.bmp");

        // Save count results
        IJ.saveAs("Results", resultsPath + "\\Results.xls");

        IJ.run("Close", "");
        outlines.close();
    }

    /*
    * Channel is the red channel that contains the cells' nuclei
    *
    * @return mask a binary mask representing the nuclei morphology for all cells
    */

    private ImagePlus generateNucleiMask(ImagePlus redChannel) {

        // Pre-processing //
        IJ.run(redChannel, "Enhance Contrast...", "saturated=0.9 equalize");

        // Processing //
        // Subtract background
        IJ.run(redChannel, "Subtract Background...", "rolling=50");

        // Threshold with black background
        IJ.setAutoThreshold(redChannel, "Default");
        Prefs.blackBackground = true;

        // Convert to mask
        IJ.run(redChannel, "Convert to Mask", "");

        // Apply closing operation to fill holes in cell nuclei
        IJ.run(redChannel, "Close-", "");

        // Apply watershed filter to separate nuclei clusters
        IJ.run(redChannel, "Watershed", "");

        return redChannel;
    }

    // Apply segmentation mask from blue channel to green channel to eliminate information outside the cells nuclei
    private ImagePlus applyMaskToGreenChannel(ImagePlus mask, ImagePlus greenChImg) {

        IJ.run(greenChImg, "8-bit", "");
        ImageCalculator ic = new ImageCalculator();
        return ic.run("Min create", mask, greenChImg);

    }

    /*
       * @param imp the segmented image of channel 3 with the channel 0 mask
       * @return the outline of all positive matches
    */
    public ImagePlus countGreenMatches(ImagePlus imp) {

        IJ.run(imp, "Enhance Contrast...", "saturated=0.9 equalize");
        IJ.setAutoThreshold(imp, "Default dark");
        Prefs.blackBackground = true;
        IJ.run(imp, "Convert to Mask", "");
        IJ.run(imp, "Analyze Particles...", "size=10-Infinity circularity=0.40-1.00 show=Outlines display clear record in_situ");

        // Return imp modified
        return imp;

    }

    public static void main(String[] args) {
        // Set the plugins.dir property to make the plugin appear in the Plugins menu
        Class<?> clazz;
        clazz = PCCNIF_Plugin.class;
        String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
        String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
        System.setProperty("plugins.dir", pluginsDir);

        // start ImageJ
        new ImageJ();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }
}