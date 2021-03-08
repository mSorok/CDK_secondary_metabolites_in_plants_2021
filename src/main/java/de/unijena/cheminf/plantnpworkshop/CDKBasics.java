/*
 * MIT License
 *
 * Copyright (c) 2021. Maria Sorokina, Aziz M. Yirik, Jonas Schaub, Christoph Steinbeck
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.plantnpworkshop;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class CDKBasics {
       /*
        * 10 randomly selected COCONUT molecules are loaded from the SDF,
        * 'COCONUTset-10.sdf' and their properties are printed.
        *
        * These molecule depictions are saved on users' desktops.
        *
        */

       /*
        * Path is set to the Desktop directory.
	*/

	public static String path= System.getProperty("user.home")+"\\Desktop\\";
	public static void main(String[] arguments) throws IOException, CDKException {
	    // Loading SDF file from resources.
            ClassLoader ClassLoader = DescriptorCalculation.class.getClassLoader();
	    File SDFile = new File(ClassLoader.getResource("COCONUTset-10.sdf").getFile());
	    IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	    IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(SDFile), builder, true);

	    //Constructors: DepictionGenerator, SmilesGenerator.
	    DepictionGenerator depiction = new DepictionGenerator();
	    SmilesGenerator sgUnique = new SmilesGenerator(SmiFlavor.Unique);
	    SmilesGenerator sgAbsolute = new SmilesGenerator(SmiFlavor.Absolute);
	    //Defining String variables.
	    String COCONUT_ID, name, uniqueSMILES, absoluteSMILES;

            //Iterating molecules
	    while (reader.hasNext()) {

	    	//Getting molecules as IAtomContainers
		IAtomContainer ac= (IAtomContainer)reader.next();

	  	//Getting its 'Name' and 'COCONUT_ID' properties.
		name = ac.getProperty("Name").toString();
		COCONUT_ID = ac.getProperty("COCONUT_ID").toString();

		//Depicting IAtomContainers. User can change depiction settings.
		depiction.withAtomColors().withCarbonSymbols().withSize(1000, 1000).withZoom(20).depict(ac).writeTo(path+COCONUT_ID+" "+name+".png");

		//Generating SMILES of IAtomContainer
		uniqueSMILES = sgUnique.create(ac);
		absoluteSMILES = sgAbsolute.create(ac);
			
		System.out.println("------------------------------");
		System.out.println("COCONUT ID and molecule name: "+COCONUT_ID+" "+name);
		System.out.println(".......");
		System.out.println("Unique SMILES: "+uniqueSMILES);
		System.out.println();
		System.out.println("Absolute SMILES: "+absoluteSMILES);
		System.out.println(".......");

		//Calculating the molecular weight.
		System.out.println(name+"'s molecular weight: "+AtomContainerManipulator.getMolecularWeight(ac));
		System.out.println(".......");
			
		//Total number of atoms.
		System.out.println("Number of atoms: "+ac.getAtomCount());

		//Number of implicit hydrogens.
		System.out.println("Implicit hydrogen count: "+AtomContainerManipulator.getImplicitHydrogenCount(ac));

		//Converting the implicit hydrogens to explicit hydrogens.
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);

		//Number of atoms after the conversion to explicit hydrogens.
		System.out.println("Number of atoms after hydrogen conversion: "+ac.getAtomCount());
	   }
	   System.out.println("------------------------------");
	   reader.close();
      }
}
