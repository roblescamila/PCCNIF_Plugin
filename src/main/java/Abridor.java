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

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import loci.common.Region;
import loci.formats.FormatException;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;

import java.io.IOException;

/**
 * A very simple example of using Bio-Formats in an ImageJ plugin.
 */
public class Abridor implements PlugIn {

    // Paths and file names
    private String dir = "C:\\Data\\PDI\\";
    private String resultsPath = dir + "results";
    private String name = "Exp035_C5_a_dcxcy3_flagcy2_bIIItubcy5_dapi_20X.czi";
    String id = dir + name;

    public void run(String arg) {
        //OpenDialog od = new OpenDialog("Open Image File...", arg);
        //String  = od.getDirectory();
        //String  = od.getFileName()

        //Open czi image
        try {
            ImporterOptions options = new ImporterOptions();
            options.setId(id);
            options.setAutoscale(true);
            options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
            options.setSplitChannels(true);
            options.isViewHyperstack();
            options.setStackOrder("XYCZT");
            ImagePlus[] imps = BF.openImagePlus(options);
            for (ImagePlus imp : imps) imp.show();
        }
        catch (FormatException exc) {
            IJ.error("Sorry, an error occurred: " + exc.getMessage());
        }
        catch (IOException exc) {
            IJ.error("Sorry, an error occurred: " + exc.getMessage());
        }
    }

    public static void main(String[] args) {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        Class<?> clazz;
        clazz = Abridor.class;
        String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
        String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
        System.setProperty("plugins.dir", pluginsDir);

        // start ImageJ
        new ImageJ();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }
}

