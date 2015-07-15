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
import ij.gui.Roi;
import ij.io.OpenDialog;
import ij.measure.ResultsTable;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.gui.GenericDialog;
import ij.process.ImageProcessor;
import loci.common.Region;
import loci.formats.FormatException;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import java.io.File;
import java.io.IOException;
import java.util.Vector;
/**
 * A very simple example of using Bio-Formats in an ImageJ plugin.
 */
public class PCCNIF_Plugin implements PlugIn {

    // Paths and file names
    String resultPath = "./result.xml";
    String resultsPath = "./results";
    String name;
    int thresholdMin = 200;
    int thresholdMax = 255;
    double circularityMin = 0.40;
    int nucleicSizeMin = 10;
    int nucleicSizeMax = 999;
    int proteinSizeMin = 10;
    int proteinSizeMax = 999;

    ResultsTable nucleos;
    ResultsTable verdes;



    public void run(String arg) {
        // Open czi image with Bio-Formats
        GenericDialog gd = new GenericDialog("Set parameters");
        gd.addNumericField("Protein threshold min", thresholdMin, 3);
        gd.addNumericField("Protein threshold max", thresholdMax, 3);
        gd.addNumericField("Nucleic circularity min", circularityMin, 3);
        gd.addNumericField("Nucleic size min", nucleicSizeMin, 3);
        gd.addNumericField("Nucleic size max", nucleicSizeMax, 3);
        gd.addNumericField("Protein size min", proteinSizeMin, 3);
        gd.addNumericField("Protein size max", proteinSizeMax, 3);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        thresholdMin = (int)gd.getNextNumber();
        thresholdMax = (int)gd.getNextNumber();
        circularityMin = (double)gd.getNextNumber();
        nucleicSizeMin = (int)gd.getNextNumber();
        nucleicSizeMax = (int)gd.getNextNumber();
        proteinSizeMin = (int)gd.getNextNumber();
        proteinSizeMax = (int)gd.getNextNumber();
        name = WindowManager.getCurrentImage().getTitle();
        name = name.substring(0,name.length()-6);

        // Close channels 2 and 3 since we do not use them
        IJ.selectWindow(name + " - C=2");
        IJ.run("Close");
        IJ.selectWindow(name + " - C=3");
        IJ.run("Close");
        /*IJ.selectWindow(name + " - C=0");
        IJ.run("View 100%");
        IJ.selectWindow(name + " - C=1");
        IJ.run("View 100%");*/

        // Generate a nuclei mask from the red channel

        ImagePlus redChannel = WindowManager.getImage(name + " - C=0");
        ImagePlus nucleiMask = generateNucleiMask(redChannel);

        // Apply the nuclei mask to the green channel
        ImagePlus greenChImg = WindowManager.getImage(name + " - C=1");
        ImagePlus outlines = countGreenMatches(greenChImg);



        // Close all unnecessary windows
        nucleiMask.changes= false;
        nucleiMask.show();
        //IJ.run("Close");

        greenChImg.changes= false;
        greenChImg.show();
        //IJ.run("Close");


        // Get positive matches
        //outlines.show();

        // Save image calculation //

        // Create results dir
        File newDir = new File(resultsPath);
        newDir.mkdirs();

        // Save image
        //IJ.run(outlines, "8-bit", "");
       // IJ.saveAs(outlines, "PNG", resultsPath + "\\" + name + " - result.bmp");

        // Save count results
        IJ.saveAs("Results", resultsPath + "\\Results.xls");

        //IJ.run("Close", "");
        //outlines.close();
    }

    /*
    * Channel is the red channel that contains the cells' nuclei
    * @return mask a binary mask representing the nuclei morphology for all cells
    */
    private ImagePlus generateNucleiMask(ImagePlus redChannel) {

        // Pre-processing //

        IJ.run(redChannel, "Enhance Contrast...", "saturated=0.9 equalize");

        // Processing //
        // Subtract background
        IJ.run(redChannel, "Subtract Background...", "rolling=50");

        // Threshold with black background
        IJ.setAutoThreshold(redChannel, "Default dark");
        Prefs.blackBackground = true;

        // Convert to mask
        IJ.run(redChannel, "Convert to Mask", "");

        // Apply closing operation to fill holes in cell nuclei
        IJ.run(redChannel, "Close-", "");

        // Apply watershed filter to separate nuclei clusters
        IJ.run(redChannel, "Watershed", "");

        IJ.run(redChannel, "Analyze Particles...", "size="+nucleicSizeMin+"-"+nucleicSizeMax+" circularity="+circularityMin+"-1.00 show=Outlines display clear record in_situ");
        nucleos = (ResultsTable) ResultsTable.getResultsTable().clone();
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

        IJ.run(imp, "8-bit", "");
        IJ.run(imp, "Enhance Contrast...", "saturated=0.9 equalize");
        IJ.run(imp, "Subtract Background...", "rolling=50");
        IJ.setThreshold(imp, thresholdMin, thresholdMax);
        Prefs.blackBackground = true;
        IJ.run(imp, "Convert to Mask", "");
        IJ.run("Set Measurements...", "area centroid center  shape decimal=3");
        IJ.run(imp, "Analyze Particles...", "size="+proteinSizeMin+"-"+proteinSizeMax+" circularity=0.00-1.00 show=Outlines display clear record in_situ");
        verdes = (ResultsTable) ResultsTable.getResultsTable().clone();
        // Return imp modified
        double[] nucleosX = nucleos.getColumnAsDoubles(ResultsTable.X_CENTER_OF_MASS);
        System.out.println("nucleosX: "+nucleosX.length);
        double[] nucleosY = nucleos.getColumnAsDoubles(ResultsTable.Y_CENTER_OF_MASS);
        double[] areas = nucleos.getColumnAsDoubles(ResultsTable.AREA);
        double[] CIs = nucleos.getColumnAsDoubles(ResultsTable.ROUNDNESS);
        double[] verdesX = verdes.getColumnAsDoubles(ResultsTable.X_CENTER_OF_MASS);
        System.out.println("verdesX: "+verdesX.length);
        double[] verdesY = verdes.getColumnAsDoubles(ResultsTable.Y_CENTER_OF_MASS);
        Vector<Vector<Integer>> result = new Vector<Vector<Integer>>();
        System.out.println(imp.getWidth()+" ,  "+imp.getHeight());
        for (int i=0; i<verdesX.length; i++){
            boolean found = false;
            int j=0;
            while(j<nucleosX.length && !found){
                if(Math.abs(nucleosX[j]-verdesX[i])<1500 && Math.abs(nucleosY[j]-verdesY[i])<1500){
                    double radio = areas[j]/Math.PI;
                    double dist = Math.pow(nucleosX[j]-verdesX[i], 2) + Math.pow(nucleosY[j]-verdesY[i], 2);

                    if(dist<= radio+ 2*radio*CIs[j]+CIs[j]*CIs[j]){
                        Vector<Integer> v = new Vector<Integer>();
                        v.add((int) (Math.min(2*nucleosX[j], imp.getWidth())));
                        v.add((int) (Math.min(2*nucleosY[j], imp.getHeight())));
                        result.add(v);
                        found = true;
                    }
                }
                j++;
            }
            if(!found) {
                ResultsTable.getResultsTable().deleteRow(i);
                ResultsTable.getResultsTable().updateResults();
            }
        }
        String resultXML = createXMLFromPoints(result);
        IJ.saveString(resultXML, resultPath);
        System.out.println(result);
        System.out.println(result.size());
        System.out.println("version 8");
        ResultsTable.getResultsTable().show("Results");
        imp.show();
        return imp;


    }

    String createXMLFromPoints(Vector<Vector<Integer>> v){
        String s = "";
        s+="<?xml version=\"1.0\" encoding=\"UTF-8\"?> <CellCounter_Marker_File> <Image_Properties> <Image_Filename>"+name+"</Image_Filename> </Image_Properties> ";
        s+="<Marker_Data> <Current_Type>1</Current_Type> <Marker_Type> <Type>1</Type> ";
        for(Vector point : v) {
            s+="<Marker > <MarkerX > "+point.elementAt(0)+" </MarkerX > <MarkerY > "+point.elementAt(1)+" </MarkerY > <MarkerZ > 1 </MarkerZ >  </Marker > ";
        }
        s+="</Marker_Type>";
        for(int i=2; i<=8; i++) s+=" <Marker_Type> <Type>"+i+"</Type> </Marker_Type> ";
        s+="</Marker_Data> </CellCounter_Marker_File> <!-- asdf -->";
        return s;
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
    }
}