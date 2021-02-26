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
 *
 * @author Jonas Schaub
 */
public class DescriptorCalculation {
    /**
     * COCONUT subset molecules are loaded from SD file and their ALogP, Petitjean number, Zagreb index, and Lipinski Rule of 5
     * violation values calculated and printed to console.
     *
     * @param args the command line arguments (none required)
     */
    public static void main(String[] args) throws FileNotFoundException, CDKException {

        //loading SD file from resources
        ClassLoader tmpClassLoader = DescriptorCalculation.class.getClassLoader();
        File tmpSDFile = new File(tmpClassLoader.getResource("COCONUTset-10.sdf").getFile());
        IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), tmpBuilder, true);

        //iterating molecules in file
        while(tmpReader.hasNext()) {
            IAtomContainer tmpMolecule = tmpReader.next();
            //reading attributes stored in SD, they are added as properties to the atom container automatically
            String tmpCOCONUTID = tmpMolecule.getProperty("COCONUT_ID");
            String tmpName = tmpMolecule.getProperty("Name");
            System.out.println("\n" + tmpName + " (" + tmpCOCONUTID + ")");

            //Petitjean number calculation
            PetitjeanNumberDescriptor tmpPetitjeanNumberDescriptor = new PetitjeanNumberDescriptor();
            tmpPetitjeanNumberDescriptor.initialise(tmpBuilder);
            DescriptorValue tmpPetitjeanNumberValue = tmpPetitjeanNumberDescriptor.calculate(tmpMolecule);
            DoubleResult tmpPetitjeanNumberResult = (DoubleResult) tmpPetitjeanNumberValue.getValue();
            System.out.println("\t" + tmpPetitjeanNumberValue.getNames()[0] + ": " + String.format("%,.2f", tmpPetitjeanNumberResult.doubleValue()));

            //Zagreb index calculation
            ZagrebIndexDescriptor tmpZagrebIndexDescriptor = new ZagrebIndexDescriptor();
            tmpZagrebIndexDescriptor.initialise(tmpBuilder);
            DescriptorValue tmpZagrebIndexValue = tmpZagrebIndexDescriptor.calculate(tmpMolecule);
            DoubleResult tmpZagrebIndexResult = (DoubleResult) tmpZagrebIndexValue.getValue();
            System.out.println("\t" + tmpZagrebIndexValue.getNames()[0] + ": " + String.format("%,.2f", tmpZagrebIndexResult.doubleValue()));

            //Lipinski Rule of 5 failures calculation
            //TODO: needs aromaticity check?
            RuleOfFiveDescriptor tmpRuleOfFiveDescriptor = new RuleOfFiveDescriptor();
            tmpRuleOfFiveDescriptor.initialise(tmpBuilder);
            DescriptorValue tmpRuleOfFiveValue = tmpRuleOfFiveDescriptor.calculate(tmpMolecule);
            IntegerResult tmpRuleOfFiveResult = (IntegerResult) tmpRuleOfFiveValue.getValue();
            System.out.println("\t" + tmpRuleOfFiveValue.getNames()[0] + ": " + tmpRuleOfFiveResult.intValue());

            //for ALogP calculation, Hs must be explicit and aromaticity detected
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpMolecule);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            //aromaticity model is constructed from electron donation model and cycle finder
            Aromaticity tmpAromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
            tmpAromaticity.apply(tmpMolecule);
            //calculates 3 values, ALogP - Ghose-Crippen LogKow, ALogP2, amr - molar refractivity
            ALOGPDescriptor tmpALogPDescriptor = new ALOGPDescriptor();
            tmpALogPDescriptor.initialise(tmpBuilder);
            //result type is an array of doubles
            DescriptorValue tmpALogPValue = tmpALogPDescriptor.calculate(tmpMolecule);
            String[] tmpALogPNames = tmpALogPValue.getNames();
            DoubleArrayResult tmpALogPResults = (DoubleArrayResult) tmpALogPValue.getValue();
            for (int i = 0; i < tmpALogPNames.length; i++) {
                System.out.println("\t" + tmpALogPNames[i] + ": " + String.format("%,.2f", tmpALogPResults.get(i)));
            }
        }
    }
}
