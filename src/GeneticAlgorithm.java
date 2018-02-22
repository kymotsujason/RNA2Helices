/*
 * GeneticAlgorithm.java by: The Kiet Ho, Jason Yue
 * CMPT166
 * Saturday, April 4, 2015
 *
 * copyright (c) 2015 The Kiet Ho, Jason Yue
 * All rights reserved
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

/**
 * This class reads in a RNA sequence from a file ("RNA.txt") and calculates all
 * possible helices associated with the sequence. Then the excess (duplicates,
 * mirrored...) helices are removed and the final results are written to a CT
 * file.
 *
 * WARNING: Writing 400+ helices will take 15+ seconds
 *
 * @author The Kiet Ho, Jason Yue
 */
public class GeneticAlgorithm {
    public static void main(String[] args) throws FileNotFoundException {
        String[] RNA = null;
        // Attempt to open "RNA.txt"; if not found, throw exception + message
        try {
            RNA = (new Scanner(new File("RNA.txt")).next()).split("");
        } catch (FileNotFoundException e) {
            System.out.println("File not found! Please make sure 'RNA.txt' exists!");
            System.exit(-1);
        }
        System.out.println("RNA length of " + RNA.length
                + " characters detected...");
        if (RNA.length > 30) {
            System.out.println("Warning! Large sequence detected! Process may "
                    + "take longer! Especially when writing to file!");
        }
        System.out.println("Begin timing...");
        // For benchmarking
        long startTime = System.nanoTime();
        List<String> permutationHelices = GenerateHelices(RNA);
        Random rand = new Random();
        // Randomly mix helix positions
        for (int i = 0; i < rand.nextInt((int) (permutationHelices.size()
                - Math.sqrt(permutationHelices.size())) + 1); i++) {
            Collections.shuffle(permutationHelices);
        }
        System.out.println(permutationHelices.size()
                + " possible Helices found...");
        permutationHelices = CleanUp(permutationHelices);
        List<String> RNAs = convertRNA(permutationHelices);
        System.out.println("Writing to CT file...");
        for (int i = 0; i < RNAs.size(); i++) {
            writeCT(RNAs.get(i), permutationHelices, i);
            if (i % 100 == 0) {
                System.out.println(".");
                System.out.println(i + "/" + RNAs.size());
            }
            System.out.print(".");
        }
        // For benchmarking
        long endTime = System.nanoTime();
        System.out.println(".");
        System.out.println("CT file succesfully created!");
        System.out.println("End timing...");
        System.out.println("Process took "
                + ((double) (endTime - startTime) / 1000000000) + " seconds!");
    }

    /**
     * (Yes, it's not the most efficient algorithm but it works quite well)
     *
     * This function attempts to build helices by matching the most common
     * canonical pairs to the sequence. Since a helix requires 3+ nucleotides
     * and base pairs, we can take 9 molecules from the sequence beginning and
     * check whether they molecules at index 0, 1, 2 match with molecules at
     * index 6, 7, 8 respectively; then we move the index up by 1 and check for
     * more matches. Once we go through the entire sequence, we increase the
     * nucleotide count to 4; thus, we take 10 molecules. Once we reach the
     * maximum amount of nucleotides, we reset the nucleotide count back to 3
     * but increase the base pair count to 4 and repeat the entire process
     * again. Once we reach the maximum amount of nucleotides again, we increase
     * the base pair count again. The base pairs are separated from the
     * nucleotide by a "|" and each pair and nucleotide are separated by "-".
     * Once we reach the max base pair count, we return an array list containing
     * all the possible helices.
     *
     * @param RNA
     * @return list of helices
     */
    public static List<String> GenerateHelices(String[] RNA) {
        System.out.println("Calculating all possible helices...");
        // Common base pairs
        String[] patterns = { "GC", "AU", "GU", "CG", "UA", "UG" };
        List<String> helices = new ArrayList<String>();
        String basePair = "", nucleotide = "";
        int lowerNucleo = 3, upperNucleo = lowerNucleo + 2, lowerPair = lowerNucleo - 3;
        int upperPair = upperNucleo + 3, pairCount = 3, nucleoCount = 2;
        while ((pairCount + (nucleoCount)) < RNA.length - 1) {
            // Since pairCount is an index, we need to multiply by 2 to
            // represent pairs
            while ((nucleoCount + (pairCount * 2)) < RNA.length - 1) {
                // Index boundary definition
                while ((RNA.length - upperPair) > 1) {
                    // We do not want to modify the original values
                    // (reliability)
                    int x = lowerPair, y = upperPair;
                    for (int i = lowerPair; i < lowerNucleo; i++) {
                        // Attempt to match each pattern
                        for (String pattern : patterns) {
                            if (pattern.matches(RNA[x] + RNA[y])) {
                                basePair += (RNA[x] + RNA[y]) + "-";
                            }
                        }
                        x++;
                        y--;
                    }
                    // Pair of 3 = 2 * 3 = 6 : "-" counts as a character; thus,
                    // 6+3 = 9 (>6 works)
                    if (basePair.length() > 6) {
                        for (int i = lowerNucleo; i <= upperNucleo; i++) {
                            nucleotide += RNA[i] + "-";
                        }
                        // Remove the dash at the end
                        basePair = basePair.substring(0, basePair.length() - 1);
                        nucleotide = nucleotide.substring(0, nucleotide.length() - 1);
                        helices.add(new String(basePair + "|" + nucleotide));
                        basePair = "";
                        nucleotide = "";
                    } else {
                        basePair = "";
                    }
                    if (upperPair < RNA.length - 1) {
                        // Keep the index moving
                        lowerNucleo++;
                        upperNucleo++;
                        lowerPair++;
                        upperPair++;
                    } else {
                        // Index is at boundary -> reset + exit while loop
                        lowerNucleo = pairCount;
                        upperNucleo = lowerNucleo + nucleoCount;
                        lowerPair = lowerNucleo - pairCount;
                        upperPair = upperNucleo + pairCount;
                    }
                }
                // Increase nucleotide count + reset with new values
                nucleoCount++;
                lowerNucleo = pairCount;
                upperNucleo = lowerNucleo + nucleoCount;
                lowerPair = lowerNucleo - pairCount;
                upperPair = upperNucleo + pairCount;
            }
            // Increase pair count + reset nucleotide count to 2 and reset with
            // new values
            pairCount++;
            nucleoCount = 2;
            lowerNucleo = pairCount;
            upperNucleo = lowerNucleo + nucleoCount;
            lowerPair = lowerNucleo - pairCount;
            upperPair = upperNucleo + pairCount;
        }
        return helices;
    }

    /**
     * This function cleans up the set of helices by removing duplicates,
     * removing helices with bases or nucleotides that are already contained in
     * other helices ("common" or non-unique helices). The HashSet library does
     * not allow duplicates, therefore we add all elements from the original set
     * (array list) of all helices to the HashSet set (which removes all
     * duplicates); then we clear the original array list and re-add all the
     * elements from HashSet back to the array list. Then we sort the array list
     * according to size, then alphabet. Next, we compare the smallest helix
     * (lowest index: 0) to the larger one (index 1) and check if the smaller is
     * contained in the larger. We then move the index and keep checking; if we
     * find a match, we remove it and check again at the same index (array lists
     * automatically shuffle indexes to fill gaps). Finally, we begin checking
     * for mirror helices by reversing the first 2 pairs and comparing, then
     * reversing the next 2 pairs and comparing and such. If we don't find a
     * match, we begin reversing the nucleotides in the same way. We then return
     * an array list of unique and organized helices.
     *
     * @param permutationHelices
     * @return list of unique and organized helices
     */
    public static List<String> CleanUp(List<String> permutationHelices) {
        System.out.println("Removing identical/excess helices...");
        // HashSets do not accept duplicates; thus, removing identical helices
        Set<String> set = new HashSet<String>();
        int allPossibleHelices = permutationHelices.size();
        set.addAll(permutationHelices);
        // Update array list with updated helices set
        permutationHelices.clear();
        permutationHelices.addAll(set);
        int removedDuplicate = permutationHelices.size();
        System.out.println(permutationHelices.size()
                + " helices remaining after removing "
                + (allPossibleHelices - permutationHelices.size())
                + " duplicates...");
        // Sort with length as priority before alphabetically
        Collections.sort(permutationHelices, lengthCompare);
        String delim = "[|]";
        int lowerIndex = 0, upperIndex = 1;
        int i = 0, match = 0;
        while (i < permutationHelices.size() - 1) {
            while (match == 0) {
                if (i < permutationHelices.size() - 1) {
                    // For comparing base pairs and nucleotides with each other
                    String[] basePair1 = permutationHelices.get(i).split(delim);
                    String[] basePair2 = permutationHelices.get(i + 1).split(delim);
                    // Checking for contained bases and nucleotides
                    // (non-unique/uninteresting)
                    if ((basePair2[0].contains(basePair1[0]))
                            || (basePair1[0].contains(basePair2[0]))) {
                        permutationHelices.remove(i);
                        match = 1;
                    } else if ((basePair2[1].contains(basePair1[1]))
                            || (basePair1[1].contains(basePair2[1]))) {
                        permutationHelices.remove(i);
                        match = 1;
                    }
                } else {
                    match = 1;
                }
                i++;
            }
            match = 0;
        }
        int removedPairs = permutationHelices.size();
        System.out.println(permutationHelices.size()
                + " helices remaining after removing "
                + (removedDuplicate - permutationHelices.size())
                + " shared pairs...");
        i = 0; match = 0;
        while (i < permutationHelices.size() - 1) {
            while (match == 0) {
                // For comparing base pairs and nucleotides with each other
                if (i < permutationHelices.size() - 1) {
                    String[] basePair1 = permutationHelices.get(i).split(delim);
                    String[] basePair2 = permutationHelices.get(i + 1).split(delim);
                    // Index boundary checking
                    if (match == 0) {
                        while (upperIndex < basePair1[0].length()) {
                            // Checking for single mirrored base pairs
                            if ((basePair2[0].contains(swap(basePair1[0].toCharArray(),
                                    lowerIndex, upperIndex)))
                                    || (swap(basePair1[0].toCharArray(), lowerIndex,
                                    upperIndex).contains(basePair2[0]))) {
                                permutationHelices.remove(i);
                                // Subtract 2 since we're adding 2 later (not that it
                                // matter)
                                upperIndex = basePair1[0].length() - 2;
                                match = 1;
                            }
                            // Add 2 to avoid switching non-matching pairs
                            lowerIndex += 2;
                            upperIndex += 2;
                        }
                        lowerIndex = 0;
                        upperIndex = 1;
                    }
                    if (match == 0) {
                        while (upperIndex < basePair2[0].length()) {
                            // Checking for single mirrored base pairs
                            if ((swap(basePair2[0].toCharArray(), lowerIndex,
                                    upperIndex).contains(basePair1[0]))
                                    || (basePair1[0].contains(swap(
                                    basePair2[0].toCharArray(), lowerIndex,
                                    upperIndex)))) {
                                permutationHelices.remove(i);
                                // Subtract 2 since we're adding 2 later (not that
                                // it matter)
                                upperIndex = basePair2[0].length() - 2;
                                match = 1;
                            }
                            // Add 2 to avoid switching non-matching pairs
                            lowerIndex += 2;
                            upperIndex += 2;
                        }
                        lowerIndex = 0;
                        upperIndex = 1;
                    }
                    if (match == 0) {
                        while (upperIndex < basePair1[1].length()) {
                            // Checking for single mirrored nucleotides
                            if ((basePair2[1].contains(swap(
                                    basePair1[1].toCharArray(), lowerIndex,
                                    upperIndex)))
                                    || (swap(basePair1[1].toCharArray(),
                                    lowerIndex, upperIndex).contains(basePair2[1]))) {
                                permutationHelices.remove(i);
                                // Subtract 2 since we're adding 2 later (not that
                                // it matter)
                                upperIndex = basePair1[1].length() - 2;
                                match = 1;
                            }
                            // Add 2 to avoid switching non-matching pairs
                            lowerIndex += 2;
                            upperIndex += 2;
                        }
                        lowerIndex = 0;
                        upperIndex = 1;
                    }
                    if (match == 0) {
                        while (upperIndex < basePair2[1].length()) {
                            // Checking for single mirrored nucleotides
                            if ((swap(basePair2[1].toCharArray(), lowerIndex,
                                    upperIndex).contains(basePair2[1]))
                                    || (basePair2[1].contains(swap(
                                    basePair2[1].toCharArray(), lowerIndex,
                                    upperIndex)))) {
                                permutationHelices.remove(i);
                                // Subtract 2 since we're adding 2 later (not that
                                // it matter)
                                upperIndex = basePair2[1].length() - 2;
                                match = 1;
                            }
                            // Add 2 to avoid switching non-matching pairs
                            lowerIndex += 2;
                            upperIndex += 2;
                        }
                        lowerIndex = 0;
                        upperIndex = 1;
                    }
                } else  {
                    match = 1;
                }
                i++;
            }
            match = 0;
        }
        System.out.println(permutationHelices.size()
                + " helices remaining after removing "
                + (removedPairs - permutationHelices.size())
                + " mirror pairs...");
        if (permutationHelices.size() >= 200 && permutationHelices.size() < 400) {
            System.out.println("Warning! Writing process will take 5+ seconds!");
        } else if (permutationHelices.size() >= 400) {
            System.out.println("Warning! Writing process will take 15+ seconds!");
        }
        return permutationHelices;
    }

    /**
     * Define a custom comparator to give priority for sorting by length
     */
    static Comparator<String> lengthCompare = new Comparator<String>() {
        @Override
        public int compare(String s1, String s2) {
            String delim = "[|]";
            String[] basePair1 = s1.split(delim);
            String[] basePair2 = s2.split(delim);
            return Integer.compare(basePair1[0].length(), basePair2[0].length());
        }
    };

    /**
     * This method is a simple character swap. It switches characters at index x
     * and index y with each other.
     *
     * @param original
     * @param lowerIndex
     * @param upperIndex
     * @return a new string with the reversed characters
     */
    public static String swap(char[] original, int lowerIndex, int upperIndex) {
        char tmp = original[lowerIndex];
        original[lowerIndex] = original[upperIndex];
        original[upperIndex] = tmp;
        return new String(original);
    }

    /**
     * This function cleans up the helices and then converts them to RNA
     * sequences. We split the nucleotides from the base pair using our
     * separator "|". We then remove the "-" from the nucleotides since it's no
     * longer needed. Next, we separate the base pairs by using our separator
     * "-". The first character of the pair (index 0) is located on the left of
     * the helix; the second of the pair (index 1) is located on the right of
     * the helix. Each additional base pair is nested within the first pair in
     * the same manner. Thus, we split the pairs by assigning each first
     * character (index 0) to a string and each additional first character on
     * the right; and similarly, we do the same with the second character (index
     * 1), but we add each additional second character on the left. Finally, we
     * put together the RNA sequence by creating a new string with the first
     * half, then the cleaned up nucleotides, then the second half. Then we add
     * the finalized RNA sequence to an array list and return the completed
     * array list once we iterate through all the helices.
     *
     * @param permutationHelices
     * @return array list of RNA sequences
     */
    public static List<String> convertRNA(List<String> permutationHelices) {
        String delim = "[|]";
        String delim2 = "[-]";
        List<String> RNAs = new ArrayList<String>();
        // First half = left side of helix before nucleotides | Second half =
        // right side of helix after nucleotides
        String[] tmp = null;
        String firstHalf = "", secondHalf = "";
        for (int i = 0; i < permutationHelices.size(); i++) {
            String[] RNA = permutationHelices.get(i).split(delim);
            // Clean up the nucleotides
            RNA[1] = RNA[1].replace("-", "");
            String[] molecules = RNA[0].split(delim2);
            for (int j = 0; j < molecules.length; j++) {
                tmp = molecules[j].split("");
                firstHalf += tmp[0];
                secondHalf = tmp[1] + secondHalf;
            }
            // Assemble RNA sequences and add to array list
            RNAs.add(firstHalf + RNA[1] + secondHalf);
        }
        return RNAs;
    }

    /**
     * This function checks whether a molecule has a pair and returns the index
     * (index + 1) of it. Using the exact same algorithm as the convertRNA
     * function, except that we do not need an array list of RNA sequences. We
     * also do not need to fully assemble the array list. Since the algorithm is
     * the same, we assume that matching pairs are in the same positions
     * (farthest left matches farthest right, 2nd farthest right matches 2nd
     * farthest right...). Thus, if the string index is less than the length of
     * the first half of the pairs, we know the matching pair is located at the
     * right of the RNA sequence. If not, we check if the string index is
     * greater than the length of the combination of the first half and the
     * nucleotide; this means the current string is an existing match with the
     * first half (refer to previous sentence). Thus, we return the positions
     * (+1) of the matched pairs. If none is matched, it means we have a
     * nucleotide and we return a 0 (meaning no pairs).
     *
     * @param RNAIndex
     * @param permutationHelices
     * @param RNA
     * @param index
     * @return the index (well, index + 1) of the matching pair
     */
    public static int baseNumber(int RNAIndex, List<String> permutationHelices,
                                 String RNA, int index) {
        String delim = "[|]", delim2 = "[-]", firstHalf_Nucleotide = "";
        // First half = left side of helix before nucleotides | Second half =
        // right side of helix after nucleotides
        String[] tmp = null;
        String firstHalf = "", secondHalf = "";
        String[] RNATypes = permutationHelices.get(RNAIndex).split(delim);
        // Clean up the nucleotides
        RNATypes[1] = RNATypes[1].replace("-", "");
        String[] molecules = RNATypes[0].split(delim2);
        for (int j = 0; j < molecules.length; j++) {
            tmp = molecules[j].split("");
            firstHalf += tmp[0];
            secondHalf = tmp[1] + secondHalf;
        }
        firstHalf_Nucleotide = firstHalf + RNATypes[1];
        /*
         * If index exists in the first half, return corresponding position (the
         * match) in second half Or, if index exists in the second half, return
         * corresponding position (the match) in first half If neither
         * conditions are met, return 0; we assume its a nucleotide.
         */
        if (index < firstHalf.length()) {
            return RNA.length() - index;
        } else if (index >= firstHalf_Nucleotide.length()) {
            if ((RNA.length() - index) <= secondHalf.length()) {
                return RNA.length() - index;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    }

    /**
     * This function writes the RNA sequence into CT format and outputs the data
     * into a file. This algorithm is based on the algorithm written by Denny
     * Chen Dai and Herbert H. Tsang.
     *
     * @param RNA
     * @param permutationHelices
     * @param RNAIndex
     */
    public static void writeCT(String RNA, List<String> permutationHelices,
                               int RNAIndex) {
        String data = "";
        for (int i = 0; i < RNA.length(); i++) {
            if (i == 0) {
                data += RNA.length() + " ENERGY = " + 0 + "\n";
            }
            // Base number; Base(A, C, G, U); Index-1; Index+1; Match Position;
            // Base number
            data += (i + 1)
                    + "\t" + RNA.charAt(i)
                    + "\t" + (i)
                    + "\t" + (i + 2)
                    + "\t"
                    + baseNumber(RNAIndex, permutationHelices, RNA, i)
                    + "\t" + (i + 1)
                    + "\n";
        }
        try {
            // Wipe file once before writing
            if (RNAIndex == 0) {
                new PrintWriter(new BufferedWriter(new FileWriter("RNA.ct",
                        false)));
            }
            PrintWriter out = new PrintWriter(new BufferedWriter(
                    new FileWriter("RNA.ct", true)));
            out.write(data);
            out.close();
        } catch (IOException e) {
            ;
        }
    }
}
