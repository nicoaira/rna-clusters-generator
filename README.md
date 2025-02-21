# RNA Triplet Dataset Generator

This project is a Rust-based RNA triplet dataset generator that mimics the original Python pipeline. It carries out the following steps:

1. **Anchor Generation**  
   A random RNA sequence is generated (using either a uniform or normal distribution), and its secondary structure is predicted via an external tool (RNAfold).

2. **Positive Sample Creation**  
   The anchor is modified to generate a positive sample. This is done by applying:
   - **Stem Modifications**: Paired base insertions or deletions are performed in stem regions (ensuring that changes happen in paired positions).
   - **Loop Modifications**: Hairpin, internal, bulge, and multiloop regions are modified with random insertions or deletions, obeying size and modification count constraints.
   - **Appending Events**: With a set probability, additional sequence segments (and their corresponding structures) are appended to one or both ends of the anchor.

3. **Negative Sample Creation**  
   A negative sample is generated by dinucleotide shuffling (preserving dinucleotide frequencies) of the anchor. Optional length variation can be applied, and the resulting sequence is folded to obtain a secondary structure.

4. **Dataset Output**  
   The generated RNA triplets (anchor, positive, and negative) are stored in a CSV file along with run metadata. Optionally, the dataset can be split into training and validation sets.

5. **Plot Generation (Optional)**  
   A specified number of triplets are selected and plotted. For each triplet, individual PNG images (one per sequence) are produced and then merged into a composite image using ImageMagick’s `montage` command.

## Usage Example

You can run the program with the following full set of command-line arguments:
