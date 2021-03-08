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

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

/**
 * Demo class for calculating selected descriptors using CDK.
 * All available descriptors can be found here:
 * <a href="http://cdk.github.io/cdk/latest/docs/api/org/openscience/cdk/qsar/IDescriptor.html"
 * >http://cdk.github.io/cdk/latest/docs/api/org/openscience/cdk/qsar/IDescriptor.html</a>
 * (All Known Implementing Classes)
 *
 * @author Jonas Schaub
 */
public class DescriptorCalculation {
    /**
     * COCONUT subset molecules are loaded from SD file and their ALogP, Petitjean number, Zagreb index, and Lipinski
     * Rule of 5 violation values calculated and printed to console.
     *
     * @param args the command line arguments (none required)
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException {

        //*loading SD file from resources*
        File tmpSDFile = new File("src/main/resources/COCONUTset-10.sdf");
        //a factory class to provide implementation independent ICDKObjects, needed for SDF parsing
        IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
        //skip: true -> erroneous entries will be skipped
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), tmpBuilder, true);

        //*iterating molecules in file*
        while(tmpReader.hasNext()) {
            IAtomContainer tmpMolecule = tmpReader.next();
            //reading attributes of molecules stored in SDF, they are added as properties to the atom container automatically
            String tmpCOCONUTID = tmpMolecule.getProperty("COCONUT_ID");
            String tmpName = tmpMolecule.getProperty("Name");
            System.out.println("\n" + tmpName + " (" + tmpCOCONUTID + ")");

            //*Petitjean number calculation*
            PetitjeanNumberDescriptor tmpPetitjeanNumberDescriptor = new PetitjeanNumberDescriptor();
            //initialise the descriptor with the specified chem object builder
            tmpPetitjeanNumberDescriptor.initialise(tmpBuilder);
            DescriptorValue tmpPetitjeanNumberValue = tmpPetitjeanNumberDescriptor.calculate(tmpMolecule);
            //DescriptorValue.getValue() returns an object implementing IDescriptorResult in general
            DoubleResult tmpPetitjeanNumberResult = (DoubleResult) tmpPetitjeanNumberValue.getValue();
            System.out.println("\t" + tmpPetitjeanNumberValue.getNames()[0] + ": " + String.format("%,.2f", tmpPetitjeanNumberResult.doubleValue()));

            //*Zagreb index calculation*
            ZagrebIndexDescriptor tmpZagrebIndexDescriptor = new ZagrebIndexDescriptor();
            tmpZagrebIndexDescriptor.initialise(tmpBuilder);
            DescriptorValue tmpZagrebIndexValue = tmpZagrebIndexDescriptor.calculate(tmpMolecule);
            DoubleResult tmpZagrebIndexResult = (DoubleResult) tmpZagrebIndexValue.getValue();
            System.out.println("\t" + tmpZagrebIndexValue.getNames()[0] + ": " + String.format("%,.2f", tmpZagrebIndexResult.doubleValue()));

            //*Lipinski Rule of 5 failures calculation*
            //preprocessing required: detection of aromaticity
            //aromaticity model is constructed from electron donation model and cycle finder
            Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
            //for detection of aromaticity, atom types must be set
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            tmpAromaticity.apply(tmpMolecule);
            RuleOfFiveDescriptor tmpRuleOfFiveDescriptor = new RuleOfFiveDescriptor();
            tmpRuleOfFiveDescriptor.initialise(tmpBuilder);
            DescriptorValue tmpRuleOfFiveValue = tmpRuleOfFiveDescriptor.calculate(tmpMolecule);
            IntegerResult tmpRuleOfFiveResult = (IntegerResult) tmpRuleOfFiveValue.getValue();
            System.out.println("\t" + tmpRuleOfFiveValue.getNames()[0] + ": " + tmpRuleOfFiveResult.intValue());

            //*ALogP calculation*
            //preprocessing required: Hs must be explicit and aromaticity detected
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
            //(done again because of now explicit Hs)
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            tmpAromaticity.apply(tmpMolecule);
            ALOGPDescriptor tmpALogPDescriptor = new ALOGPDescriptor();
            tmpALogPDescriptor.initialise(tmpBuilder);
            //calculates 3 values: ALogP (Ghose-Crippen LogKow), ALogP2, amr (molar refractivity)
            DescriptorValue tmpALogPValue = tmpALogPDescriptor.calculate(tmpMolecule);
            String[] tmpALogPNames = tmpALogPValue.getNames();
            //result type is an array of 3 doubles
            DoubleArrayResult tmpALogPResults = (DoubleArrayResult) tmpALogPValue.getValue();
            //iterating the 3 calculated values
            for (int i = 0; i < tmpALogPNames.length; i++) {
                System.out.println("\t" + tmpALogPNames[i] + ": " + String.format("%,.2f", tmpALogPResults.get(i)));
            }
        } //end of while()
    } //end of main()
} //end of class
