package com.example.dyrka2
//chatgpt do listy kodonów

/**
 * Abstract class for other nucleotide sequence classes. Contains a String identifier, stores data
 * as a String, has a Set of allowed Chars and an Int length. Returns error for invalid characters.
 */
open class Sequence(
    val identifier: String,
    var data: String,
    private val validChars: Set<Char>
){
    val length: Int
        get() = data.length

    init{
        require(data.all {it in validChars}){"Invalid characters in sequence."}
    }

    /**
     * Functions returns the sequence as FASTA type string.
     * Input:
     * -
     * Return:
     * sequence data as FASTA (String)
     */
    override fun toString(): String{
        return ">$identifier\n$data"
    }

    /**
     * Function performs a given mutation on the sequence in a chosen spot.
     * Input:
     * -position (Int) - position in the sequence for the mutation to occur
     * -value (Char) - mutation that happens in the chosen spot
     * Returns:
     * -doesn't return anything, just changes the data String
     * Exceptions:
     * -if invalid character (not in the VALID_CHARS set)
     * -if position is out of range
     */
    open fun mutate(position: Int, value: Char) {
        require(value in validChars){"Invalid character for mutation."}
        require(position in data.indices){"Mutation is out of range."}

        val charArray = data.toCharArray()
        charArray[position] = value
        data = String(charArray)
    }

    /**
     * Function finds the chosen String motif in the sequence.
     * Input:
     * -motif (String) - motif to search for
     * Returns:
     * -motif index (Int) - index of the motif, -1 if not found
     */
    fun findMotif(motif: String): Int{
        return data.indexOf(motif)
    }
}

/**
 * Class that inherits from the abstract class Sequence
 */
class DNASequence(identifier: String, data: String) : Sequence(identifier, data, setOf('A', 'T', 'G', 'C')){
    /**
     * Function returns the complementary strand for the DNA sequence.
     * Input:
     * -
     * Returns:
     * result (String) - complementary strand
     */
    fun complement(): String{
        var result = ""
        for (char in data) {
            result += when(char) {
                'A' -> 'T'
                'T' -> 'A'
                'G' -> 'C'
                'C' -> 'G'
                else -> char
            }
        }
        return result
    }

    /**
     * Function transcribes DNA to RNA.
     * Input:
     * -
     * Returns:
     * - RNA sequence (RNASequence) - sequence RNA for the given DNA strand.
     */
    fun transcribe(): RNASequence{
        val transcribed = data.replace('T', 'U')
        return RNASequence(identifier, transcribed)
    }
}

/**
 * Class that inherits from the abstract class Sequence
 */
class RNASequence(identifier: String, data: String) : Sequence(identifier, data, setOf('A', 'U', 'G', 'C')){
    /**
     * Function transcribes the sequence to create a protein sequence.
     * Input:
     * -
     * Returns:
     * - protein sequence (ProteinSequence) - protein sequence object after transcribing the RNA
     */
    fun transcribe(): ProteinSequence{
        var protein = ""
        for (i in 0 until data.length step 3){
            if (i + 2 < data.length){
                val c1 = data[i]
                val c2 = data[i + 1]
                val c3 = data[i + 2]

                val codon = "" + c1 + c2 + c3
                //poprosiłem chatgpt o wypisanie tej konwersji, bo nie ma mowy, że będę to robił ręcznie
                val aminoAcid = when (codon) {
                    "GCU", "GCC", "GCA", "GCG" -> 'A'
                    "UGU", "UGC" -> 'C'
                    "GAU", "GAC" -> 'D'
                    "GAA", "GAG" -> 'E'
                    "UUU", "UUC" -> 'F'
                    "GGU", "GGC", "GGA", "GGG" -> 'G'
                    "CAU", "CAC" -> 'H'
                    "AUU", "AUC", "AUA" -> 'I'
                    "AAA", "AAG" -> 'K'
                    "UUA", "UUG", "CUU", "CUC", "CUA", "CUG" -> 'L'
                    "AUG" -> 'M'
                    "AAU", "AAC" -> 'N'
                    "CCU", "CCC", "CCA", "CCG" -> 'P'
                    "CAA", "CAG" -> 'Q'
                    "CGU", "CGC", "CGA", "CGG", "AGA", "AGG" -> 'R'
                    "UCU", "UCC", "UCA", "UCG", "AGU", "AGC" -> 'S'
                    "ACU", "ACC", "ACA", "ACG" -> 'T'
                    "GUU", "GUC", "GUA", "GUG" -> 'V'
                    "UGG" -> 'W'
                    "UAU", "UAC" -> 'Y'
                    "UAA", "UAG", "UGA" -> '|'  //stop
                    else -> 'X' //błąd
                }
                protein += aminoAcid
            }
        }
        return ProteinSequence(identifier, protein)
    }
}

/**
 * Class that inherits from the abstract class Sequence
 */
class ProteinSequence(identifier: String, sequence: String) : Sequence(identifier, sequence, ('A'..'Z').toSet())

fun main(){
    var s1=DNASequence("ABC", "ATCGAAACGGCTATC")
    println("Your DNA sequence is: ${s1.data}")
    check(s1.data=="ATCGAAACGGCTATC"){"Wrong sequence data!"}

    println("Your DNA identifier is: ${s1.identifier}")
    check(s1.identifier=="ABC"){"Wrong identifier!"}

    println("FASTA format for your sequence:\n$s1")
    check(s1.toString()==">ABC\nATCGAAACGGCTATC"){"Bad FASTA format!"}

    println("Length of the sequence is: ${s1.length}")
    check(s1.length==15){"Wrong length!"}

    s1.mutate(2, 'A')
    println("Your mutated sequence is:\n$s1")
    check(s1.data=="ATAGAAACGGCTATC"){"Mutation is wrong!"}

    val complement = s1.complement()
    println("Complement of the sequence: $complement")
    check(complement == "TATCTTTGCCGATAG") { "Complement is incorrect!" }

    val rna = s1.transcribe()
    println("Transcribed RNA sequence:\n$rna")
    check(rna.data == "AUAGAAACGGCUAUC") { "Transcription to RNA is incorrect!" }

    println("Position of your motif ACG is: ${s1.findMotif("ACG")}")
    check(s1.findMotif("ACG")==6){"Wrong spot for the motif!"}
    println("Position of your motif TTT is: ${s1.findMotif("TTT")}")
    check(s1.findMotif("TTT")==-1){"The motif doesn't exist in this sequence!"}

    val protein = rna.transcribe()
    println("Translated protein sequence:\n$protein")
    check(protein.data == "IETAI") {"Transcription to protein is incorrect!"}

    try{
        var s2 = DNASequence("DEF", "AUCG")
    }catch (e: IllegalArgumentException){
        println("\n${e.message}")
    }
    try{
        var s2 = DNASequence("DEF", "ATCG")
        s2.mutate(0, 'X')
    }catch (e: IllegalArgumentException){
        println("\n${e.message}")
    }
    try{
        var s2 = DNASequence("DEF", "ATCG")
        s2.mutate(10, 'T')
    }catch (e: IllegalArgumentException){
        println("\n${e.message}")
    }
}