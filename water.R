#------------------------------------------------------------------------------
# Create a data frame that resembles the format of a .hb2 file, see definitions
# in Table 1 of https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html
# Thankfully, this is possible because .hb2 files are highly structured. Note
# that we need to run HBPLUS on a pdb file. From now on, unless specified, we
# will call this data frame 'hbonds'
#
# Input:
#   pdb.name: four-letter pdb name, no extension required, e.g. '1AD2'
# 
# Return a data frame that includes numeric and string variables. From now on,
# unless specified, we will call this data frame 'hbonds'
#------------------------------------------------------------------------------
process_hb2 <- function(pdb.name) {
  # Read the respective .hb2 file line by line
  # con stands for contents
  hb2.file <- file(paste(pdb.name, '.hb2', sep = ''), "r")
  con <- c()
  while ( TRUE ) {
    line = readLines(hb2.file, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    con <- c(con, line)
  }
  close(hb2.file)
  con <- con[9:length(con)]    # skip comments from line 1 to 8
  
  # An example of an entry in .hb2 file, read top down for the column number in 
  # Table 1 of https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html
  # 
  # 00000000011111111112222222222333333333344444444445555555555666666666677777777778
  # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
  # A0069-VAL N   A2001-HOH O   3.27 MH  -2 -1.00 163.8  2.28  -1.0  -1.0     1
  #
  # donor           :  A0069 (whitespace removed)
  # donor_amino     :  VAL   (whitespace removed)
  # donor_type      :  N     (whitespace removed)
  # acceptor        :  A2001 (whitespace removed)
  # acceptor_amino  :  HOH   (whitespace removed)
  # acceptor_type   :  O     (whitespace removed)
  # da_dist         :  3.27
  # cat_da          :  MH    (whitespace removed)
  # ... so on, refer to Table 1 ...
  # 
  n <- length(con)
  hbonds <- data.frame(donor = character(n), donor_amino = character(n), 
                       donor_type = character(n), acceptor = character(n), 
                       acceptor_amino = character(n), acceptor_type = character(n), 
                       da_dist = numeric(n), cat_da = character(n), 
                       num_aas = numeric(n), dist = numeric(n), 
                       dha_angle = numeric(n), ha_dist = numeric(n), 
                       h_a_aa_angle = numeric(n), d_a_aa_angle = numeric(n), 
                       stringsAsFactors = FALSE)
  for (i in 1:length(con)) {
    ro <- substring(con[i], first = c(1, 7, 10, 15, 21, 24, 28, 34, 37, 41, 47, 53, 59, 65), 
                    last  = c(5, 9, 13, 19, 23, 27, 32, 35, 39, 45, 51, 57, 63, 69))
    ro[2] <- trimws(ro[2])
    ro[3] <- trimws(ro[3])
    ro[5] <- trimws(ro[5])
    ro[6] <- trimws(ro[6])
    ro[7] <- as.numeric(ro[7])
    ro[9] <- as.numeric(ro[9])
    ro[10] <- as.numeric(ro[10])
    ro[11] <- as.numeric(ro[11])
    ro[12] <- as.numeric(ro[12])
    ro[13] <- as.numeric(ro[13])
    ro[14] <- as.numeric(ro[14])
    hbonds[i, ] <- ro 
  }
  return(hbonds)
}

#------------------------------------------------------------------------------
# Check if a donor or acceptor amino is a nucleotide or not as defined in 
# Table 2 from https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html
# 
# Input
#   x: amino acid code
# 
# Return TRUE of FALSE
#------------------------------------------------------------------------------
is.nucleotide <- function(x) {
  return(x %in% c('C', 'A', 'U', 'G', 'T', 'ATP'))
}

#------------------------------------------------------------------------------
# Check if a donor or acceptor amino is an amino acid or not as defined in 
# Table 2 from https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html
#
# Input
#   x: amino acid code
# 
# Return TRUE of FALSE
#------------------------------------------------------------------------------
is.std.amino.acid <- function(x) {
  std.amino.acids <- c('ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX',
                       'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                       'THR', 'TRP', 'TYR', 'VAL')
  return(x %in% std.amino.acids)
}

#------------------------------------------------------------------------------
# Process each hydrogen bond of hbonds data frame, and indicate if the hydrogen
# bond is between a nucleotide and an amino acid or not. To obtain how many of
# such hydrogen bonds, try sum(nucle.aacid(hbonds)).
#
# Input
#   hbonds: a data frame produced by function 'process_hb2'
# 
# Return a vector of TRUE/FALSE, where TRUE indicates a hydrogen bond between a
# nucleotide and an amino acid; and FALSE otherwise.
#------------------------------------------------------------------------------
nucle.aacid <- function(hbonds) {
  case1 <- is.nucleotide(hbonds$donor_amino)     & is.std.amino.acid(hbonds$acceptor_amino)
  case2 <- is.std.amino.acid(hbonds$donor_amino) & is.nucleotide(hbonds$acceptor_amino)
  return(case1 | case2)
}

#------------------------------------------------------------------------------
# Create a network base on a hbonds data frame, where each vertex is the chain
# ID together with residue number, e.g. A0069, called ID for simplicity. From 
# each vertex, we can obtain a list of adjacent vertices that may be waters, 
# nucleotides, or amino acids. This is a graph data structure.
#
# Input
#   hbonds: a data frame produced by function 'process_hb2'
# 
# Return a vector in which each element is a list that represent a vertex.
# Each vertex will contains information about which adjacent vertices are waters,
# nucleotides, and/or amino acids 
#------------------------------------------------------------------------------
hb.network <- function(hbonds) {
  # Initialize the vertices
  vertices <- unique(c(hbonds$donor, hbonds$acceptor))    # obtain all IDs
  vertices <- sort(vertices)    # I like to sort to debug easier, but still works fine without
  n.vertices <- length(vertices)
  network <- vector("list", n.vertices)   # create a vector of list
  names(network) <- vertices
  
  # Process each vertex and obtain adjacent vertices
  for (i in 1:n.vertices) {
    # Filter all hydrogen bonds  
    id <- vertices[i]
    df <- hbonds[hbonds$donor == id | hbonds$acceptor == id, ]
    
    hoh <- unique(c(df$donor[df$donor_amino == 'HOH'], df$acceptor[df$acceptor_amino == 'HOH']))
    hoh <- hoh[hoh != id]
    nucle <- unique(c(df$donor[is.nucleotide(df$donor_amino)], df$acceptor[is.nucleotide(df$acceptor_amino)]))
    nucle <- nucle[nucle != id]
    aacid <- unique(c(df$donor[is.std.amino.acid(df$donor_amino)], df$acceptor[is.std.amino.acid(df$acceptor_amino)]))
    aacid <- aacid[aacid != id]
    network[[i]]$id <- id
    network[[i]]$hoh <- hoh
    network[[i]]$nucle <- nucle
    network[[i]]$aacid <- aacid
  }
  
  return(network)
}

#------------------------------------------------------------------------------
# Identify single water bridges base on a hbonds data frame. A single bridge is
# formed when a water molecule mediates a nucleotide and an amino acid like this
#
#         nucleotide ----- water ----- amino acid
# 
# the reverse is the same. The water forms at least one hydrogen bond with the
# nucleotide, and also forms at least one hydrogen bond with the amino acid. The
# number of hydrogen bonds between them is not taken into account. For instace, if
# there are 2 hydrogen bonds formed between the water and the nucleotide, and 3 
# hydrogen bonds formed between the water and the amino acid; collectively, that 
# counts as 1 single water bridge. However, those hydrogen bonds are called
# water-mediated hydrogen bonds; and we will count them seperately at a different
# time.
#
# An output of this function will look like this
#
#            nucle   hoh aacid
#            R0128 L1003 L0056
#            R0129 L1001 H0106
#            R0129 L1002 L0049
#            R0129 L1001 L0049
#
# Note that although the pairs R0129-H0106 and R0129-L0049 are mediated by the same
# water L1001, they form 2 distinct single water bridges. So each row in 
# 
# Input
#   hbonds    : a data frame produced by function 'process_hb2'
# 
# Return a data frame where each row corresponds to a single water bridge.
#------------------------------------------------------------------------------
find.swb <- function(hbonds) {
  # Obtain all the nucleotides, and filter out those that are possible to form 
  # a single water brigde. To identify those nucleotides, I process each nucleotide,
  # and obtain all hydrogen bonds that involves such nucleotide. If these hydrogen
  # bonds do not associate with any water, then there is no way this nucleotide can
  # form a single water bridge (how can I form a water bridge if I am not in touch with
  # any water???)
  # The function may still work just fine without the filtering part, but probably slower.
  # Filtering should be fast, as it used 'sapply' function, which works exactly as a loop,
  # but much faster than a loop.
  nucles <- unique(c(hbonds$donor[is.nucleotide(hbonds$donor_amino)], hbonds$acceptor[is.nucleotide(hbonds$acceptor_amino)]))
  nucles <- sort(nucles)    # I like sorting to debug easier
  nucles <- nucles[!sapply(nucles, function(nucle) {
    df <- hbonds[hbonds$donor == nucle | hbonds$acceptor == nucle, ]
    all(df$donor_amino != 'HOH') & all(df$acceptor_amino != 'HOH')
  })]
  
  # Do the same thing as what I did for the nucleotides, but with amino acids now
  aacids <- unique(c(hbonds$donor[is.std.amino.acid(hbonds$donor_amino)], hbonds$acceptor[is.std.amino.acid(hbonds$acceptor_amino)]))
  aacids <- sort(aacids)    # I like sorting to debug easier
  aacids <- aacids[!sapply(aacids, function(aacid) {
    df <- hbonds[hbonds$donor == aacid | hbonds$acceptor == aacid, ]
    all(df$donor_amino != 'HOH') & all(df$acceptor_amino != 'HOH')
  })]
  
  network <- hb.network(hbonds)    # get the network
  waters <- matrix(, nrow = 0, ncol = 3)   # an empty matrix to store output, will convert to data frame later
  
  # Basically, here, I start with nucleotides (the reverse is the same), and for each nucleotide, I check each 
  # amino acid if I can go from the nucleotide to a water to the amino acid, each exactly in one step.
  for (nucle in nucles) {    
    for (aacid in aacids) {
      # recall that I sort vertices in 'hb.network'. List in R does not work well with character indices, so 
      # I have to convert to numerical indices
      ref.nucle <- which(names(network) == nucle)    
      for (hoh in network[[ref.nucle]]$hoh) {    # network[[ref.nucle]]$hoh : adjacent waters of the nucleotide
        ref.hoh <- which(names(network) == hoh)
        if (aacid %in% network[[ref.hoh]]$aacid) {    # network[[ref.hoh]]$aacid : adjacent amino acids of the water
          waters <- rbind(waters, c(nucle, hoh, aacid))
        }
      }
    }
  }
  
  # just some formatting to make the output informative
  waters <- as.data.frame(waters)
  names(waters) <- c('nucle', 'hoh', 'aacid')
  return(waters)
}

#------------------------------------------------------------------------------
# Identify double water bridges base on a hbonds data frame. A double bridge is
# formed when 2 connected water molecules mediate a nucleotide and an amino acid 
# like this
#
#         nucleotide ----- water 1 ----- water 2 ----- amino acid
# 
# the reverse is the same. Water 1 forms at least one hydrogen bond with the
# nucleotide, and also forms at least one hydrogen bond with water 2. Water 2 
# forms at least one hydrogen bonds with the amino acid. The number of hydrogen
# bonds between them is not taken into account. For instace, if
# there are 2 hydrogen bonds formed between wate r1 and the nucleotide, and 3 
# hydrogen bonds formed between water 2 and the amino acid; collectively, that 
# counts as 1 double water bridge. I don't think we came up with a terminology
# for hydrogen bonds that involves in a double water bridge yet, I will leave it
# here for you (sorry!)
#
# An output of this function will look like this
#
#        nucle  hoh1  hoh2 aacid
#        R0128 R1049 L1003 L0056
#        R0129 R1125 H1001 H0110
#        R0129 R1208 H1001 H0110
#        R0129 R1186 L1047 L0049
#        R0129 L1005 L1002 L0049
#
# Note that seems like the pair R0129-H0110 is repeated, but the intermediate waters
# can be different, so R1125 and H1001 form a double water bridge between R0129 and H0110,
# 
# Input
#   hbonds    : a data frame produced by function 'process_hb2'
# 
# Return a data frame where each row corresponds to a double water bridge.
#------------------------------------------------------------------------------
find.dwb <- function(hbonds) {
  # as specified in 'find.swb'
  nucles <- unique(c(hbonds$donor[is.nucleotide(hbonds$donor_amino)], hbonds$acceptor[is.nucleotide(hbonds$acceptor_amino)]))
  nucles <- sort(nucles)
  nucles <- nucles[!sapply(nucles, function(nucle) {
    df <- hbonds[hbonds$donor == nucle | hbonds$acceptor == nucle, ]
    all(df$donor_amino != 'HOH') & all(df$acceptor_amino != 'HOH')
  })]
  
  # as specified in 'find.swb'
  aacids <- unique(c(hbonds$donor[is.std.amino.acid(hbonds$donor_amino)], hbonds$acceptor[is.std.amino.acid(hbonds$acceptor_amino)]))
  aacids <- sort(aacids)
  aacids <- aacids[!sapply(aacids, function(aacid) {
    df <- hbonds[hbonds$donor == aacid | hbonds$acceptor == aacid, ]
    all(df$donor_amino != 'HOH') & all(df$acceptor_amino != 'HOH')
  })]
  
  network <- hb.network(hbonds)
  waters <- matrix(nrow = 0, ncol = 4)
  for (nucle in nucles) {
    for (aacid in aacids)  {
      ref.nucle <- which(names(network) == nucle)
      for (hoh1 in network[[ref.nucle]]$hoh) {
        ref.hoh1 <- which(names(network) == hoh1)
        for (hoh2 in network[[ref.hoh1]]$hoh) {
          ref.hoh2 <- which(names(network) == hoh2) 
          if (aacid %in% network[[ref.hoh2]]$aacid) {
            waters <- rbind(waters, c(nucle, hoh1, hoh2, aacid))
          }
        }
      }
    }
  }
  
  waters <- as.data.frame(waters)
  names(waters) <- c('nucle', 'hoh1', 'hoh2', 'aacid')
  return(waters)
}

#------------------------------------------------------------------------------
# Count or return the unique waters base on a hbonds data frame. These are
# chain IDs together with residue numbers, called IDs in short, e.g. R1245, R1149.
# Note that the function does not return the hydrogen bonds associated with these
# unique waters. That would be done seperately.
#
# Input
#   hbonds: a data frame produced by function 'process_hb2'
#   count: TRUE/FALSE, whether the count or the list of unique waters should be
#          returned. TRUE as default.
# 
# Return the number of unique waters if count = TRUE, or the list of unique waters
# if count = FALSE
#------------------------------------------------------------------------------
unique.waters <- function(hbonds, count = TRUE) {
  all.waters <- c(hbonds$donor[hbonds$donor_amino == 'HOH'], hbonds$acceptor[hbonds$acceptor_amino == 'HOH'])
  all.waters <- unique(all.waters)
  if (!count) {
    return(all.waters)
  }
  return(length(all.waters))
}

#------------------------------------------------------------------------------
# Count or return the unique amino acids base on a hbonds data frame. These are
# chain IDs together with residue numbers, called IDs in short.
# Note that the function does not return the hydrogen bonds associated with these
# unique amino acids. That would be done seperately.
#
# Input
#   hbonds: a data frame produced by function 'process_hb2'
#   count: TRUE/FALSE, whether the count or the list of unique amino acids should be
#          returned. TRUE as default.
# 
# Return the number of unique amino acids if count = TRUE, or the list of unique 
# amino acids if count = FALSE
#------------------------------------------------------------------------------
unique.amino.acid <- function(hbonds, count = TRUE) {
  all.proteins <- c(hbonds$donor[is.std.amino.acid(hbonds$donor_amino)], hbonds$acceptor[is.std.amino.acid(hbonds$acceptor_amino)])
  all.proteins <- unique(all.proteins)
  if (!count) {
    return(all.proteins)
  }
  return(length(all.proteins))
}

#------------------------------------------------------------------------------
# Count or return the unique nucleotides base on a hbonds data frame. These are
# chain IDs together with residue numbers, called IDs in short.
# Note that the function does not return the hydrogen bonds associated with these
# unique nucleotides. That would be done seperately.
#
# Input
#   hbonds: a data frame produced by function 'process_hb2'
#   count: TRUE/FALSE, whether the count or the list of unique nucleotides should be
#          returned. TRUE as default.
# 
# Return the number of unique nucleotides if count = TRUE, or the list of unique 
# nucletides if count = FALSE
#------------------------------------------------------------------------------
unique.nucleotide <- function(hbonds, count = TRUE) {
  all.rnas <- c(hbonds$donor[is.nucleotide(hbonds$donor_amino)], hbonds$acceptor[is.nucleotide(hbonds$acceptor_amino)])
  all.rnas <- unique(all.rnas)
  if (!count) {
    return(all.rnas)
  }
  return(length(all.rnas))
}

#------------------------------------------------------------------------------
# Create a data frame that resembles the interfacing residue table from PDBePisa.
# When lauching PDBePisa, put the pdb name, click Analyze, click View Details,
# the interfacing residue table will locate at the very bottom. A dated version 
# of the table can also be found at https://www.ebi.ac.uk/pdbe/docs/Tutorials/workshop_tutorials/PDBepisa.pdf
# on the top of page 6. The function will need the .xml file of the table. Although
# we can name that .xml anything, to be consistent, I like to name it using PDB name.
# An output of this function will look like this
#
#   id structure hsdc        asa        bsa bsa_score solv_energy
#    1      GLU8       99.622000   0.000000         0  0.00000000
#    1      ARG9       75.412900   0.000000         0  0.00000000
#    1     PRO10       44.692200   0.000000         0  0.00000000
#
# At the time of writing this function, only 'structure', 'bsa', and 'bsa_score'
# are of interest, but a thorough processing of this table is done anyway for 
# future use. A bsa_score > 0 indicates that the corresponding structure is in the
# interfacial region. Note that the chain IDs (A and B) used by PDBePisa are not
# consistent with the chain IDs used by .hb2 file of HBPLUS. This could be confusing,
# but I managed to map that back to the hbonds data frame (.hb2 file) in function
# 'interfacial.hbond'
#
# Input
#   pdb.name: four-letter pdb name, no extension required, e.g. '1AD2'
#
# Return a data frame that represent the interfacing residue table.
#------------------------------------------------------------------------------
extract.residues <- function(pdb.name) {
  pdb.name <- paste(pdb.name, '.xml', sep = '')
  require(xml2)    # we need this to process a .xml file
  
  # An XML document has a tree structure. The one from PDBePisa has 2 big nodes,
  # each has a number of child nodes. Our goal is to visit each node and extract
  # the information. The contents of an XML file are usually strings, so for
  # number information, I want to make sure that I convert to numeric in the data frame.
  xml.doc <- read_xml(pdb.name, options = '')
  res <- NULL
  for (j in 1:2) {    # 2 big nodes
    node <- xml_children(xml.doc)[j]
    residues <- xml_children(node)
    for (i in 1:length(residues)) {    # children of each big node
      df <- data.frame(id = NA, structure = character(1), hsdc = character(1), asa = numeric(1), 
                       bsa = numeric(1), bsa_score = numeric(1), solv_energy = numeric(1))
      df$id <- j
      residue <- xml_children(residues[i])
      
      # Here, I want to normalize the format of the structure column in [amino code][residue number with 1 to 4 digits]
      structure.node <- residue[1]
      x <- trimws(xml_text(structure.node))
      x <- unlist(strsplit(x, split = '[A-Z]:'))[2]
      x <- gsub('\\s+', '', x)
      df$structure <- x
      
      hsdc.node <- residue[2]
      x <- trimws(xml_text(hsdc.node))
      df$hsdc <- x
      
      asa.node <- residue[3]
      x <- xml_text(asa.node)
      x <- as.numeric(x)
      df$asa <- x
      
      bsa.node <- residue[4]
      x <- xml_text(bsa.node)
      x <- as.numeric(x)
      df$bsa <- x
      
      bsa.score.node <- residue[5]
      x <- xml_text(bsa.score.node)
      x <- as.numeric(x)
      df$bsa_score <- x
      
      solv.energy.node <- residue[6]
      x <- xml_text(solv.energy.node)
      x <- as.numeric(x)
      df$solv_energy <- x
      
      res <- rbind(res, df)
    }
  }
  return(res)
}

#------------------------------------------------------------------------------
# Indicate in the hbonds data frame which hydrogen bonds are associated with the
# interfacial regions determined by the yellow region in PDBePisa (structures with
# bsa_score > 0). Due to the inconsistency of chain ID in the .xml file and .hb2 file
# donor and acceptor in the hbonds data frame are formatted into 
# [amino code][residue number with 1 to 4 digits], e.g. ARG9, PRO123
#
# Input
#   pdb.name: four-letter pdb name, no extension required, e.g. '1AD2'
#   regions: a vector of structures that reside in the interfacial region.
#
# Return the indices of hydrogen bonds associated with the interfacial regions.
#------------------------------------------------------------------------------
interfacial.hbond <- function(hbonds, regions) {
  n <- nrow(hbonds)
  if (n == 0) {
    return(integer(0))
  }
  
  df <- data.frame(don = character(n), acc = character(n))
  for (i in 1:n) {
    b <- hbonds$donor[i]
    b <- unlist(strsplit(b, split = '[A-Z]'))[2]
    b <- as.numeric(b)
    a <- hbonds$donor_amino[i]
    df$don[i] <- paste(a, b, sep = '')
    
    b <- hbonds$acceptor[i]
    b <- unlist(strsplit(b, split = '[A-Z]'))[2]
    b <- as.numeric(b)
    a <- hbonds$acceptor_amino[i]
    df$acc[i] <- paste(a, b, sep = '')
  }
  return(which(df$don %in% regions | df$acc %in% regions))
}
