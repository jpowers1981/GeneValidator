
##
# Codon bias for one amino acid
class CodonBias

  attr_accessor :aa_short_code
  attr_accessor :aa_code
  attr_accessor :aa_name
  attr_accessor :codons # array of possible codons
  attr_accessor :bias # counts the no of apparitions of each codon

  def initialize(aa_short_code, aa_code, aa_name, codons)
    @aa_short_code = aa_short_code.upcase
    @aa_code = aa_code.upcase
    @aa_name = aa_name.capitalize
    @codons = codons.map{|codon| codon.upcase}
    @bias = Hash.new(0)
  end

  def add_amino_acid(codon)
    if codons.include?(codon)
      bias[codon.upcase] = bias[codon.upcase] + 1
    end
  end

  def add_codon_bias(codon_bias)
    codon_bias.bias.each do |codon, freq|
      if codons.include?(codon)
        @bias[codon.upcase] = @bias[codon.upcase] + freq
      end
    end
  end

  def total_count
    total_count = 0
    codons.each do |codon|
      total_count = total_count + bias[codon]
    end
    return total_count
  end

  def get_percentage(codon)
    unless codons.include?(codon)
      return nil
    end
    return bias[codon] / (total_count + 0.0)
  end

  def print_bias
    print "#{@aa_name}(#{@aa_short_code}): "
    codons.each do |codon|
      print "#{codon}:#{get_percentage(codon).round(2)} "
    end
    puts ""
  end

end


class AllCodonsBias < Array

  attr_reader :ala, :asx, :cys, :asp, :glu, :phe, :gle, :his, :ile, :lys, :leu, :met
  attr_reader :asl, :pro, :gln, :arg, :ser, :thr, :val, :trp, :tyr, :glx 
  attr_reader :map_name_codon
 
  def initialize

    # compute the codon bias for the prediction
    @ala = CodonBias.new("A", "ALA", "Alanine", ["GCA", "GCC", "GCG", "GCT"])
    @asx = CodonBias.new("B", "ASX","Asparagine", ["AAC", "AAT", "GAC", "GAT"])
    @cys = CodonBias.new("C", "CYS", "Cysteine", ["TGC", "TGT"])
    @asp = CodonBias.new("D", "ASP", "Aspartic acid", ["GAC", "GAT"])
    @glu = CodonBias.new("E", "GLU", "Glutamic acid", ["GAA", "GAG"])
    @phe = CodonBias.new("F", "PHE", "Phenylalanine", ["TTC", "TTT"])
    @gle = CodonBias.new("G", "GLE", "Glycine", ["GGA", "GGC", "GGG", "GGT"])
    @his = CodonBias.new("H", "HIS", "Histidine", ["CAC", "CAT"])
    @ile = CodonBias.new("I", "ILE", "Isoleucine", ["ATA", "ATC", "ATT"])
    @lys = CodonBias.new("K", "LYS", "Lysine", ["AAA", "AAG"])
    @leu = CodonBias.new("L", "LEU", "Leucine", ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"])
    @met = CodonBias.new("M", "MET", "Methionine", ["ATG"])
    @asl = CodonBias.new("N", "ASL", "Asparagine", ["AAC", "AAT"])
    @pro = CodonBias.new("P", "PRO", "Proline", ["CCA", "CCC", "CCG", "CCT"])
    @gln = CodonBias.new("Q", "GLN", "Glutamine", ["CAA", "CAG"])
    @arg = CodonBias.new("R", "ARG", "Arginine", ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"])
    @ser = CodonBias.new("S", "SER", "Serine", ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"])
    @thr = CodonBias.new("T", "THR", "Threonine", ["ACA", "ACC", "ACG", "ACT"])
    @val = CodonBias.new("V", "VAL", "Valine", ["GTA", "GTC", "GTG", "GTT"])
    @trp = CodonBias.new("W", "TRP", "Tryptophan", ["TGG"])
    @tyr = CodonBias.new("Y", "TYR", "Tyrosine", ["TAC", "TAT"])
    @glx = CodonBias.new("Z", "GLX", "Glutamine", ["CAA", "CAG", "GAA", "GAG"])
    @nn = CodonBias.new("NN", "NN", "Unknowm", [])

    @map_name_codon = {}
    @map_name_codon["A"] = @ala
    @map_name_codon["B"] = @asx
    @map_name_codon["C"] = @cys
    @map_name_codon["D"] = @asp
    @map_name_codon["E"] = @glu
    @map_name_codon["F"] = @phe
    @map_name_codon["G"] = @gle
    @map_name_codon["H"] = @his
    @map_name_codon["I"] = @ile
    @map_name_codon["K"] = @lys
    @map_name_codon["L"] = @leu
    @map_name_codon["M"] = @met
    @map_name_codon["N"] = @asl
    @map_name_codon["P"] = @pro
    @map_name_codon["Q"] = @gln
    @map_name_codon["R"] = @arg
    @map_name_codon["S"] = @ser
    @map_name_codon["T"] = @thr
    @map_name_codon["V"] = @val
    @map_name_codon["W"] = @trp
    @map_name_codon["Y"] = @tyr
    @map_name_codon["Z"] = @glx
    @map_name_codon["NN"] = @nn

    push @ala
    push @asx
    push @cys
    push @asp
    push @glu
    push @phe
    push @gle
    push @his
    push @ile
    push @lys
    push @leu
    push @met
    push @asl
    push @pro
    push @gln
    push @arg
    push @ser
    push @thr
    push @val
    push @trp
    push @tyr
    push @glx
  end

  def get_map_codon_aa
    map_codon_aa = {"GCA"=>@ala,"GCC"=>@ala,"GCG"=>@ala, "GCT"=>@ala,
                    "AAC"=>@asx, "AAT"=>@asx, "GAC"=>@asx, "GAT"=>@asx,
                    "TGC"=>@cys, "TGT"=>@cys,
                    "GAC"=>@asp, "GAT"=>@asp,
                    "GAA"=>@glu, "GAG"=>@glu,
                    "TTC"=>@phe, "TTT"=>@phe,
                    "GGA"=>@gle, "GGC"=>@gle, "GGG"=>@gle, "GGT"=>@gle,
                    "CAC"=>@his, "CAT"=>@his,
                    "ATA"=>@ile, "ATC"=>@ile, "ATT"=>@ile,
                    "AAA"=>@lys, "AAG"=>@lys,
                    "CTA"=>@leu, "CTC"=>@leu, "CTG"=>@leu, "CTT"=>@leu, "TTA"=>@leu, "TTG"=>@leu,
                    "ATG"=>@met,
                    "AAC"=>@asl, "AAT"=>@asl,
                    "CCA"=>@pro, "CCC"=>@pro, "CCG"=>pro, "CCT"=>pro,
                    "CAA"=>@gln, "CAG"=>@gln,
                    "AGA"=>@arg, "AGG"=>@arg, "CGA"=>@arg, "CGC"=>@arg, "CGG"=>@arg, "CGT"=>@arg,
                    "AGC"=>@ser, "AGT"=>@ser, "TCA"=>@ser, "TCC"=>@ser, "TCG"=>@ser, "TCT"=>@ser,
                    "ACA"=>@thr, "ACC"=>@thr, "ACG"=>@thr, "ACT"=>@thr,
                    "GTA"=>@val, "GTC"=>@val, "GTG"=>@val, "GTT"=>@val,
                    "TGG"=>@trp,
                    "TAC"=>@tyr, "TAT"=>@tyr,
                    "CAA"=>@glx, "CAG"=>@glx, "GAA"=>@glx, "GAG"=>@glx}

    return map_codon_aa
  end
  
  # input: AllCodonsBias
  def update_codon_bias(codon_bias)
    codon_bias.each do |codon| 
      @map_name_codon[codon.aa_short_code].add_codon_bias(codon)
    end
  end

end

