package org.enzymeanalyzer;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.io.GenbankReaderHelper;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import java.io.*;
import java.util.*;

public class Main {
    public static void main(String[] args) throws Exception {

        String filePath = args[0].replace("\\", "/");

        File dnaFile = new File(filePath);

        LinkedHashMap<String, DNASequence> dnaSequences;
        dnaSequences = GenbankReaderHelper.readGenbankDNASequence(dnaFile);
        DNASequence sequence = dnaSequences.get(args[1]);

        BufferedReader reader = new BufferedReader(new FileReader(args[2].replace("\\", "/")));
        ArrayList<String> allEnzymes = new ArrayList<>();
        String enzLine;
        while ((enzLine = reader.readLine()) != null) {
            allEnzymes.add(enzLine);
        }
        reader.close();

        for (String enz : allEnzymes) {
            ArrayList<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> featureArray = new ArrayList<>();
            for (FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> i : sequence.getFeaturesByType("CDS")) {
                try {
                    if (i.getQualifiers().get("EC_number").get(0).getValue().equals(enz)) {
                        featureArray.add(i);
                    }
                } catch (NullPointerException ignored) {
                    
                }
            }

            if (featureArray.size() == 0) {
                throw new Exception("No such enzyme in the genome");
            }

            PrintWriter out = new PrintWriter(args[3].replace("\\", "/") + "EnzymeInfo_" + enz + ".txt");
            out.println("- - - Amino acids statistics - - -");
            for (FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound> seq : featureArray) {
                ProteinSequence protSeq = new ProteinSequence(seq.getQualifiers().get("translation").get(0).getValue());
                AminoAcidCompoundSet allAminoAcids = new AminoAcidCompoundSet();
                out.println("protein id: " + seq.getQualifiers().get("protein_id").get(0).getValue().replaceAll("[.](,?.*)", ""));
                out.print("Polar: ");
                out.print(protSeq.countCompounds(allAminoAcids.getCompoundForString("S")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("T")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("C")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("N")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("Q")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("Y"))
                );
                out.print(", ");
                out.print("Hydrophobic: ");
                out.println(protSeq.countCompounds(allAminoAcids.getCompoundForString("G")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("A")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("V")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("L")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("I")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("P")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("F")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("M")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("W"))
                );
                out.print("Aliphatic: ");
                out.print(protSeq.countCompounds(allAminoAcids.getCompoundForString("G")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("A")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("V")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("L")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("I"))
                );
                out.print(", ");
                out.print("Aromatic: ");
                out.print(protSeq.countCompounds(allAminoAcids.getCompoundForString("F")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("Y")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("W"))
                );
                out.print(", ");
                out.print("Heterocyclic: ");
                out.println(protSeq.countCompounds(allAminoAcids.getCompoundForString("P")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("H")) +
                        protSeq.countCompounds(allAminoAcids.getCompoundForString("W"))
                );
                out.println();
            }
            HashMap<String, String> enzymePathwaysMap = EnzymePathways.get(enz);
            out.println("- - - Enzyme pathways - - -");
            for (Map.Entry<String, String> entry : enzymePathwaysMap.entrySet()) {
                String key = entry.getKey();
                String value = entry.getValue();
                out.println("Pathway: " + key + "  Description: " + value);
            }
            out.close();
        }
    }
}