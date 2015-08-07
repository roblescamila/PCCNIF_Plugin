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
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import java.io.File;
import java.util.Vector;
/**
 * A very simple example of using Bio-Formats in an ImageJ plugin.
 */
public class PCCNIF_Plugin implements PlugIn {

    // Paths and file names
    String resultPath = "./result.xml";
    String resultsPath = "./results";
    String name;
    double enhancecontrastsaturated = 0.4;
    int substractbackgroundrolling = 50;
    int nucleithresholdMin = 50;
    int nucleithresholdMax = 255;
    double circularityMin = 0.40;
    int nucleicSizeMin = 40;
    int nucleicSizeMax = 4000;

    int proteinthresholdMin = 100;
    int proteinthresholdMax = 255;
    int proteinSizeMin = 40;
    int proteinSizeMax = 4000;
    int nucleichannel = 1;
    int proteinchannel = 4;
    int proteinthreshold = 100;
    boolean processnuclei = false;
    boolean processprotein = false;
    boolean proteincounting = false; // false método de rectángulo, true método de conteo
    ImageProcessor nucleiimage;
    ImageProcessor proteinimage;

    ResultsTable nucleos;
    ResultsTable verdes;




    public void run(String arg) {
        // Open czi image with Bio-Formats
        GenericDialog gd = new GenericDialog("Set parameters");

        gd.addNumericField("Ehance contrast saturation", enhancecontrastsaturated, 3);
        gd.addNumericField("Substract Background rolling", substractbackgroundrolling, 3);
        gd.addNumericField("Nuclei threshold min", nucleithresholdMin, 3);
        gd.addNumericField("Nuclei threshold max", nucleithresholdMax, 3);
        gd.addNumericField("Nucleic circularity min", circularityMin, 3);
        gd.addNumericField("Nucleic size min", nucleicSizeMin, 3);
        gd.addNumericField("Nucleic size max", nucleicSizeMax, 3);

        gd.addNumericField("Protein threshold min", proteinthresholdMin, 3);
        gd.addNumericField("Protein threshold max", proteinthresholdMax, 3);
        gd.addNumericField("Protein size min", proteinSizeMin, 3);
        gd.addNumericField("Protein size max", proteinSizeMax, 3);

        gd.addNumericField("Nuclei channel", nucleichannel, 1);
        gd.addNumericField("Protein channel", proteinchannel, 1);
        gd.addCheckbox("Preprocess Nuclei", processnuclei);
        gd.addCheckbox("Preprocess Protein", processprotein);
        gd.addCheckbox("Protein counting",proteincounting);
        //    gd.add
        gd.showDialog();
        if (gd.wasCanceled()) return;
        enhancecontrastsaturated = (double)gd.getNextNumber();
        substractbackgroundrolling = (int)gd.getNextNumber();
        nucleithresholdMin = (int)gd.getNextNumber();
        nucleithresholdMax = (int)gd.getNextNumber();
        circularityMin = (double)gd.getNextNumber();
        nucleicSizeMin = (int)gd.getNextNumber();
        nucleicSizeMax = (int)gd.getNextNumber();

        proteinthresholdMin = (int)gd.getNextNumber();
        proteinthresholdMax = (int)gd.getNextNumber();
        proteinSizeMin = (int)gd.getNextNumber();
        proteinSizeMax = (int)gd.getNextNumber();
        nucleichannel = (int)gd.getNextNumber();
        proteinchannel = (int)gd.getNextNumber();
        processnuclei = gd.getNextBoolean();
        processprotein = gd.getNextBoolean();
        proteincounting = gd.getNextBoolean();

        // Generate a nuclei mask from the red channel
        name = WindowManager.getCurrentImage().getTitle();
        ImageStack is =  WindowManager.getCurrentImage().getStack();
        nucleiimage = is.getProcessor(nucleichannel);
        proteinimage = is.getProcessor(proteinchannel);

        ImagePlus nucleiip = new ImagePlus();
        nucleiip.setProcessor(nucleiimage);
        nucleiip.show();

        ImagePlus proteinip = new ImagePlus();
        proteinip.setProcessor(proteinimage);
        proteinip.show();
        if (processnuclei)
        {
            nucleiip = generateNucleiMask(nucleiip);
        }

        IJ.run("Set Measurements...", "area centroid center bounding shape decimal=3");
        String APparameters = "size=" + nucleicSizeMin + "-" + nucleicSizeMax + " circularity=" + circularityMin + "-1.00 show=Outlines display clear record in_situ";
        IJ.run(nucleiip, "Analyze Particles...", APparameters);
        nucleos = (ResultsTable) ResultsTable.getResultsTable().clone();

        if (processprotein)
        {
            proteinip = processproteinChannel(proteinip);
        }
        if (proteincounting)
        {
            ImagePlus outlines = countGreenMatches(proteinip);
            outlines.show();
        }
        else
            contarpositivos(proteinimage);
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

        IJ.run(redChannel, "8-bit", "");
        // Pre-processing //
        IJ.run(redChannel, "Enhance Contrast...", "saturated=" + enhancecontrastsaturated + " equalize");

        // Processing //
        // Subtract background
        IJ.run(redChannel, "Subtract Background...", "rolling="+substractbackgroundrolling);

        // Threshold with black background
        IJ.setAutoThreshold(redChannel, "Default dark");
        IJ.setThreshold(redChannel,nucleithresholdMin,nucleithresholdMax);
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

    public ImagePlus processproteinChannel(ImagePlus imp) {
        IJ.run(imp, "8-bit", "");
        IJ.run(imp, "Enhance Contrast...", "saturated=0.9 equalize");
        IJ.run(imp, "Subtract Background...", "rolling=50");
        IJ.setThreshold(imp, proteinthresholdMin, proteinthresholdMax);
        Prefs.blackBackground = true;
        IJ.run(imp, "Convert to Mask", "");
        return imp;
    };

    public int contarpositivos(ImageProcessor imp) {
        int cantidadpositivos = 0;
        for (int i=0;i<ResultsTable.getResultsTable().size()-1;i++)
        {
            // obtener el ROI de los núcleos
            int area = (int) ResultsTable.getResultsTable().getValue("Area",i);
            int x = (int) ResultsTable.getResultsTable().getValue("BX",i);
            int y = (int) ResultsTable.getResultsTable().getValue("BY",i);
            int width = (int) ResultsTable.getResultsTable().getValue("Width",i);
            int height = (int) ResultsTable.getResultsTable().getValue("Height",i);
            // get all type-specific pixels (fast)
            // in this example, a ByteProcessor
            int cantidadmayorathreshold = 0;
            for (int k = x; k < x+width; k++)
                for (int l = y; l < y+height; l++) {
                    // Java has no unsigned 8-bit data type, so we need to perform Boolean arithmetics
                    int value = imp.getPixel(k,l);
                    if (value > proteinthreshold)
                        cantidadmayorathreshold++;
                }
            if (cantidadmayorathreshold/area > 0.5)
                cantidadpositivos++;
        }
        return cantidadpositivos;
    }

    public ImagePlus countGreenMatches(ImagePlus imp) {

        IJ.run("Set Measurements...", "area centroid center shape decimal=3");
        IJ.run(imp, "Analyze Particles...", "size="+proteinSizeMin+"-"+proteinSizeMax+" circularity=0.00-1.00 show=Outlines display clear record in_situ");
        verdes = (ResultsTable) ResultsTable.getResultsTable().clone();
        // Return imp modified
        double[] nucleosX = nucleos.getColumnAsDoubles(ResultsTable.X_CENTER_OF_MASS);
        double[] nucleosY = nucleos.getColumnAsDoubles(ResultsTable.Y_CENTER_OF_MASS);
        double[] areas = nucleos.getColumnAsDoubles(ResultsTable.AREA);
        double[] CIs = nucleos.getColumnAsDoubles(ResultsTable.ROUNDNESS);
        double[] verdesX = verdes.getColumnAsDoubles(ResultsTable.X_CENTER_OF_MASS);
        double[] verdesY = verdes.getColumnAsDoubles(ResultsTable.Y_CENTER_OF_MASS);
        Vector<Vector<Integer>> result = new Vector<Vector<Integer>>();
        for (int i=0; i<verdesX.length; i++){
            boolean found = false;
            int j=0;
            while(j<nucleosX.length && !found){
                if(Math.abs(nucleosX[j]-verdesX[i])<2*nucleicSizeMax && Math.abs(nucleosY[j]-verdesY[i])<2*nucleicSizeMax){
                    double radio = areas[j]/Math.PI;
                    double dist = Math.pow(nucleosX[j]-verdesX[i], 2) + Math.pow(nucleosY[j]-verdesY[i], 2);

                    if(dist<= radio+ 2*radio*CIs[j]+CIs[j]*CIs[j]){
                        Vector<Integer> v = new Vector<Integer>();
                        v.add((int) nucleosX[j]);
                        v.add((int) nucleosY[j]);
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
        ResultsTable.getResultsTable().show("Results");
        imp.show();
        return imp;


    }

    String createXMLFromPoints(Vector<Vector<Integer>> v){
        String s = "";
        s+="<?xml version=\"1.0\" encoding=\"UTF-8\"?> <CellCounter_Marker_File> <Image_Properties> <Image_Filename>Counter Window - "+name+"</Image_Filename> </Image_Properties> ";
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