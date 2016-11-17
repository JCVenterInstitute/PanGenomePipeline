# Exception raised when a parse error happens in processing a Newick tree
class NewickParseError < RuntimeError
end

# Represents a token (substring with meaning) in a Newick parse
class NewickToken
  # semantic meaning of token (label, weight, etc.)
  attr_reader :type
  # string value of token
  attr_reader :value

  def initialize(type, value)
    @type = type
    @value = value
  end

end


# Splits a Newick tree string into tokens that NewickTree uses
class NewickTokenizer

  def initialize(str)
    @str = str
    @pos = 0
  end

  # returns the next character in the string and updates position
  def nextChar
    if (@pos < @str.size)
      c = @str[@pos].chr
      @pos += 1
      return c
    else
      return nil
    end
  end

  # returns the next token in the string and updates position
  def nextToken
    c = nextChar
    if (c == " " || c == "\n" || c == "\r")
      return nextToken
    elsif (c == "(" || c == ")" || c == ',')
      return NewickToken.new("SYMBOL", c)
    elsif (c == ":")
      if (@str.index((/([0-9|\.|\-|e|E]+)/), @pos) == @pos)
        @pos += $1.length
        return NewickToken.new("WEIGHT", $1)
      else
        raise NewickParseError, "Illegal weight at pos #{@pos} of #{@str}"
      end
    elsif (c == "'")
      if (@str.index(/(\'[^\']*\')/, @pos - 1) == @pos - 1)
        @pos += $1.length - 1
        return NewickToken.new("LABEL", $1)
      else
        raise NewickParseError, "Illegal label at pos #{@pos} of #{@str}"
      end
    elsif (@str.index(/([^,():]+)/, @pos - 1) == @pos - 1)
      @pos += $1.length - 1
      return NewickToken.new("LABEL", $1)
    end
  end

  # returms the next token in the string without changing position
  def peekToken
    origPos = @pos
    token = nextToken
    @pos = origPos
    return token
  end

end

# Represents a single node in a NewickTree
class NewickNode
  # parent node of node
  attr :parent, true
  # edge length of node
  attr :edgeLen, true
  # name of node
  attr :name, true
  # child nodes of node
  attr_reader :children
  # x position of node
  attr :x, true
  # y position of node
  attr :y, true

  def initialize(name, edgeLen)
    @parent = nil
    @name = name
    @edgeLen = edgeLen
    @children = []
  end

  # adds child node to list of children and sets child's parent to self
  def addChild(child)
    child.parent = self
    @children.push(child)
  end

  # removes child node from list of children and sets child's parent to nil
  def removeChild(child)
    @children.delete(child)
    child.parent = nil
  end

  # returns string representation of node
  def to_s(showLen = true, bootStrap = "node")
    s = ""
    if (!leaf?)
      s += "("
      @children.each do |child|
        s += child.to_s(showLen, bootStrap)
        s += "," if (child != @children.last)
      end
      s += ")"
    end
    if (leaf? || bootStrap == "node")
      s += @name
    end
    s += ":#{@edgeLen}" if (showLen && @edgeLen != 0)
    if (!leaf? && name.to_i > 0 && bootStrap == "branch")
      s += ":#{name}"
    end
    return s
  end

  # returns array of names of leaves (taxa) that are contained in the node
  def taxa(bootstrap = false)
    taxa = []
    if (!leaf?)
      taxa.push(@name) if (bootstrap)
      @children.each do |child|
        child.taxa.each do |taxon|
          taxa.push(taxon)
        end
      end
    else
      taxa.push(name)
    end
    return taxa.sort
  end

  # returns array of leaves (taxa) are contained in the node
  def leaves
    nodes = []
    descendants.each do |node|
      nodes.push(node) if (node.leaf?)
    end
    return nodes
  end

  # returns array of non leaves (taxa) that are contained in the node
  def intNodes
    nodes = []
    descendants.each do  |child|
      nodes.push(child) if (!child.leaf?)
    end
    return nodes
  end

  # returns node with given name, or nil if not found
  def findNode(name, exact=false)
    found = nil
    if (exact && @name==name) || (!exact && @name =~/#{name}/)
      found = self
    else
      @children.each do |child|
        found = child.findNode(name, exact)
        break if found
      end
    end
    return found
  end

  # reverses the parent-child relationship (used in rerooting tree)
  def reverseChildParent
    return if (@parent.nil?)
    oldParent = @parent
    oldParent.removeChild(self)
    if (!oldParent.parent.nil?)
      oldParent.reverseChildParent
    end
    addChild(oldParent)
    oldParent.edgeLen = @edgeLen
    @edgeLen = 0
  end

  # True if given node is child (or grandchild, etc.) of self. False otherwise
  def include?(node)
    while(node.parent != nil)
      return true if (node.parent == self)
      node = node.parent
    end
    return false
  end

  # True if node has no children (and therefore is a leaf)
  def leaf?
    if (@children.empty?)
      return true
    else
      return false
    end
  end

  # returns array of all descendant nodes
  def descendants
    descendants = []
    @children.each do |child|
      descendants.push(child)
      child.descendants.each do |grandchild|
        descendants.push(grandchild)
      end
    end
    return descendants
  end

  # return array of all sibling nodes
  def siblings
    siblings = []
    if (parent.nil?)
      return siblings
    else
      @parent.children.each do |child|
        siblings.push(child) if (child!=self)
      end
      return siblings
    end
  end

  # reorders descendant nodes alphabetically and by size
  def reorder
    return if (@children.empty?)
    @children.sort! {|x, y| x.name <=> y.name}
    @children.each do |child|
      child.reorder
    end
    return self
  end



  # returns the last common ancestor node of self and given node
  def lca(node)
    if (self.include?(node))
      return self
    elsif (node.include?(self))
      return node
    else
      return @parent.lca(node)
    end
  end

  # returns the distance to the ancestor node
  def distToAncestor(ancestor)
    dist = 0
    node = self
    while(node != ancestor)
      dist += node.edgeLen
      node = node.parent
    end
    return dist
  end


  # returns number of nodes to the ancestor node
  def nodesToAncestor(ancestor)
    if (!ancestor.include?(self))
      return nil
    elsif (ancestor == self)
      return 0
    elsif (ancestor == @parent)
      return 1
    else
      return 1 + @parent.nodesToAncestor(ancestor)
    end
  end


  # returns number of nodes to other node
  def nodesToNode(node)
    lca = lca(node)
    if (lca == self)
      return node.nodesToAncestor(self)
    elsif (lca == node)
      return nodesToAncestor(node)
    else
      return nodesToAncestor(lca) + node.nodesToAncestor(lca)
    end
  end

  # calculates node Y positions
  def calcYPos
    ySum = 0
    @children.each do |child|
      ySum += child.y
    end
    @y = ySum / @children.size
  end

  # calculates node X positions
  def calcXPos
    if (parent.nil?)
      @x = 0
    else
      #@edgeLen = 1 if (@edgeLen == 0)
      @x = parent.x + @edgeLen
    end
    if (!leaf?)
      @children.each do |child|
        child.calcXPos
      end
    end
  end

  # returns the maximum X value in node
  def xMax
    xMax = 0
    children.each do |child|
      xMax = child.edgeLen if (child.x > xMax)
    end
    return xMax
  end

  # returns the maximum Y value in node
  def yMax
    yMax = 0
    children.each do |child|
      yMax = child.y if (child.y > yMax)
    end
    return yMax
  end

  # returns the minimum Y value in node
  def yMin
    yMin = 1e6
    children.each do |child|
      yMin = child.y if (child.y < yMin)
    end
    return yMin
  end

end

class NewickTree
  attr_reader :root
  def initialize(treeString)
    tokenizer = NewickTokenizer.new(treeString)
    @root = buildTree(nil, tokenizer)
  end

  # create new NewickTree from tree stored in file
  def NewickTree.fromFile(fileName)
    treeString = ""
    inFile = File.new(fileName)
    inFile.each do |line|
      treeString += line.chomp
    end
    inFile.close
    treeString.gsub!(/\[[^\]]*\]/,"") # remove comments before parsing
    return NewickTree.new(treeString)
  end

  # internal function used for building tree structure from string
  def buildTree(parent, tokenizer)
    while (!(token = tokenizer.nextToken).nil?)
      if (token.type == "LABEL")
        name = token.value
        edgeLen = 0
        if (tokenizer.peekToken.type == "WEIGHT")
          edgeLen = tokenizer.nextToken.value.to_f
        end
        node = NewickNode.new(name, edgeLen)
        return node
      elsif (token.value == "(")
        node = NewickNode.new("", 0)
        forever = true
        while (forever)
          child = buildTree(node, tokenizer)
          node.addChild(child)
          break if tokenizer.peekToken.value != ","
          tokenizer.nextToken
        end
        if (tokenizer.nextToken.value != ")")
          raise NewickParseError, "Expected ')' but found: #{token.value}"
        else
          peek = tokenizer.peekToken
          if (peek.value == ")" || peek.value == "," || peek.value == ";")
            return node
          elsif (peek.type == "WEIGHT")
            node.edgeLen = tokenizer.nextToken.value.to_f
            return node
          elsif (peek.type == "LABEL")
            token = tokenizer.nextToken
            node.name = token.value
            if (tokenizer.peekToken.type == "WEIGHT")
              node.edgeLen = tokenizer.nextToken.value.to_f
            end
            return node
          end
        end 
      else
        raise NewickParseError, 
        "Expected '(' or label but found: #{token.value}"
      end
    end
  end

  # return string representation of tree
  def to_s(showLen = true, bootStrap = "node")
    return @root.to_s(showLen, bootStrap) + ";"
  end

  # write string representation of tree to file
  def write(fileName, showLen = true, bootStrap = "node")
    file = File.new(fileName, "w")
    file.print @root.to_s(showLen, bootStrap) + ";\n"
    file.close
  end

  # reorders leaves alphabetically and size
  def reorder
    @root.reorder
    return self
  end

  # renames nodes and creates an alias file, returning aliased tree and hash
  def alias(aliasFile = nil, longAlias = false)
    ali = Hash.new
    aliF = File.new(aliasFile, "w") if (!aliasFile.nil?)
    if (longAlias)
      taxon = "SEQ" + "0"* taxa.sort {|x,y| x.length <=> y.length}.last.length
    else
      taxon =  "SEQ0000001"
    end
    @root.descendants.each do |node|
      if (node.name != "" && node.name.to_i == 0)
        ali[taxon] = node.name
        aliF.printf("%s\t%s\n", taxon, node.name) if (!aliasFile.nil?)
        node.name = taxon.dup
        taxon.succ!
      end
    end
    aliF.close if (!aliasFile.nil?)
    return self, ali
  end

  # renames nodes according to alias hash
  def unAlias(aliasNames)
    @root.descendants.each do |node|
      node.name = aliasNames[node.name] if (!aliasNames[node.name].nil?)
    end
    return self
  end

  # renames nodes according to inverse alias hash 
  def reAlias(aliasNames)
    @root.descendants.each do |node|
      aliasNames.keys.each do |key|
        node.name = key if (aliasNames[key] == node.name)
      end
    end
    return self
  end

  # return array of all taxa in tree
  def taxa
    return @root.taxa
  end

  # returns a 2D hash of pairwise distances on tree
  def distanceMatrix(verbose = nil)
    dMatrix = Hash.new
    leaves = root.leaves
    num = (leaves.size)**2
    numpercent = num/100
    count = 0
    
    leaves.each do |leaf1|
      dMatrix[leaf1.name] = Hash.new
      leaves.each do |leaf2|
        count += 1
        if (verbose && count % numpercent == 0)
          STDERR.printf("Processing %s %s (%d of %d)\n", leaf1.name, leaf2.name, count, num)
        end 
        if (leaf1 == leaf2)
          dMatrix[leaf1.name][leaf2.name] = 0.0
        else
          lca = leaf1.lca(leaf2)
          dMatrix[leaf1.name][leaf2.name] = leaf1.distToAncestor(lca) + leaf2.distToAncestor(lca)
        end
      end
    end
    return dMatrix
  end

  # returns lists of clades different between two trees
  def compare(tree)
    tree1 = self.dup.unroot
    tree2 = tree.dup.unroot

    diff1 = []
    diff2 = []
    if (tree1.taxa == tree2.taxa)
      clades1 = tree1.clades
      clades2 = tree2.clades
      clades1.each do |clade|
        if (!clades2.include?(clade))
          diff1.push(clade)
        end
      end
      clades2.each do |clade|
        if (!clades1.include?(clade))
          diff2.push(clade)
        end
      end
    else
      raise NewickParseError, "The trees have different taxa!"
    end
    return diff1, diff2
  end

  # return node with the given name, exactly if second argument true
  def findNode(name, exact=false)
    return @root.findNode(name, exact)
  end

  # unroot the tree
  def unroot
    if (@root.children.size != 2)
      return self # already unrooted
    end
    left, right = @root.children
    left, right = right, left if (right.leaf?) # don't uproot leaf side   
    left.edgeLen += right.edgeLen
    right.children.each do |child|
      @root.addChild(child)
    end
    @root.removeChild(right)
    return self
  end

  # root the tree on a given node
  def reroot(node)
    unroot
    left = node
    right = left.parent
    right.removeChild(node)
    right.reverseChildParent
    if (left.edgeLen != 0)
      right.edgeLen = left.edgeLen / 2.0
      left.edgeLen = right.edgeLen
    end
    @root = NewickNode.new("", 0)
    @root.addChild(left)
    @root.addChild(right)
    return self
  end

  # returns the two most distant leaves and their distance apart
  def mostDistantLeaves
    greatestDist = 0
    dist = Hash.new
    org1, org2 = nil, nil
    @root.leaves.each do |node1|
      @root.leaves.each do |node2|
        dist[node1] = Hash.new if dist[node1].nil?
        dist[node2] = Hash.new if dist[node2].nil?
        next if (!dist[node1][node2].nil?)
        lca = node1.lca(node2)
        dist[node1][node2] = node1.distToAncestor(lca) + 
        node2.distToAncestor(lca)
        dist[node2][node1] = dist[node1][node2]
        if (dist[node1][node2] > greatestDist)
          org1 = node1
          org2 = node2
          greatestDist = dist[node1][node2]
        end
      end
    end
    return org1, org2, greatestDist
  end

  # add EC numbers from alignment
  def addECnums(alignFile)
    ec = Hash.new
    File.new(alignFile).each do |line|
      if (line =~ /^>/)
        definition = line.chomp[1..line.length]
        name = definition.split(" ").first
        if (definition =~ /\[EC:([0-9|\.]*)/)
          ec[name] = name + "_" + $1
        end
      end
    end
    unAlias(ec)
  end

  # root the tree on midpoint distance
  def midpointRoot
    unroot
    org1, org2, dist = mostDistantLeaves
    midDist = dist / 2.0
    return self if (midDist == 0)
    if (org1.distToAncestor(@root) > org2.distToAncestor(@root))
      node = org1
    else
      node = org2
    end
    distTraveled = 0
    while(!node.nil?)
      distTraveled += node.edgeLen
      break if (distTraveled >= midDist)
      node = node.parent
    end
    oldDist = node.edgeLen
    left, right = node, node.parent
    right.removeChild(node)
    right.reverseChildParent
    left.edgeLen = distTraveled - midDist
    right.edgeLen = oldDist - left.edgeLen
    @root = NewickNode.new("", 0)
    @root.addChild(left)
    @root.addChild(right)
    return self
  end

  # returns array of arrays representing the tree clades
  def clades(bootstrap = false)
    clades = []
    @root.descendants.each do |clade|
      clades.push(clade.taxa(bootstrap)) if (!clade.children.empty?)
    end
    return clades
  end

  # add bootstrap values (given in clade arrays) to a tree
  def addBootStrap(bootClades)
    @root.descendants.each do |clade|
      next if clade.leaf?
      bootClades.each do |bClade|
        boot, rest = bClade.first, bClade[1..bClade.size - 1]
        if (rest == clade.taxa ) # same clade found
          clade.name = boot
        end
      end
    end
  end

  # return array of arrays of taxa representing relatives at each level
  def relatives(taxon)
    node = findNode(taxon)
    if (node.nil?)
      return nil
    else
      relatives = []
      while(!node.parent.nil?)
        relatives.push(node.parent.taxa - node.taxa)
        node = node.parent
      end
      return relatives
    end
  end


  # Fixes PHYLIP's mistake of using branch lengths and not node values
  def fixPhylip
    @root.descendants.each do |child|
      br = child.edgeLen.to_i
      child.edgeLen = 0
      if (br > 0 && !child.leaf?)
        child.name = br.to_s
      end
    end
  end


  # calculates leaf node positions (backwards from leaves, given spacing)
  def calcPos(yUnit)
    yPos = 0.25
    @root.reorder
    leaves = @root.leaves.sort {|x, y| x.nodesToNode(y) <=> y.nodesToNode(x)}
    leaves.each do |leaf|
      leaf.y = yPos
      yPos += yUnit
    end
    nodes =  @root.intNodes.sort{|x, y| y.nodesToAncestor(@root) <=> 
      x.nodesToAncestor(@root)}
      nodes.each do |node|
        node.calcYPos
      end
      @root.calcYPos
      @root.calcXPos
      nodes =  @root.intNodes.sort{|x, y| x.nodesToAncestor(@root) <=> 
        y.nodesToAncestor(@root)}
        nodes.each do |node|
          @root.calcXPos # (forwards from root)
        end
      end

      # function to generate gi link to ncbi for draw, below
      def giLink(entry)
        ncbiLink = "http://www.ncbi.nlm.nih.gov/entrez/"
        protLink = "viewer.fcgi?db=protein&val="
        if (entry =~ /^gi[\_]*([0-9]*)/ || entry =~ /(^[A-Z|0-9]*)\|/)
          return ncbiLink + protLink + $1
        else
          return nil
        end
      end

      # returns PDF representation of branching structure of tree
      def draw(pdfFile, boot="width", linker = :giLink, labelName = false,
        highlights = Hash.new, brackets = nil, rawNames = false, goodBoot = 75)
        colors = {"goldenrod4"=>[139, 105, 20], "lightcyan"=>[224, 255, 255], "rosybrown1"=>[255, 193, 193], "sienna1"=>[255, 130, 71], "lavender"=>[230, 230, 250], "rosybrown2"=>[238, 180, 180], "wheat1"=>[255, 231, 186], "linen"=>[250, 240, 230], "beige"=>[245, 245, 220], "sienna2"=>[238, 121, 66], "paleturquoise1"=>[187, 255, 255], "pink1"=>[255, 181, 197], "wheat2"=>[238, 216, 174], "violetred"=>[208, 32, 144], "palegreen"=>[152, 251, 152], "gray50"=>[127, 127, 127], "rosybrown3"=>[205, 155, 155], "sienna3"=>[205, 104, 57], "paleturquoise2"=>[174, 238, 238], "green1"=>[0, 255, 0], "pink2"=>[238, 169, 184], "wheat3"=>[205, 186, 150], "lightskyblue1"=>[176, 226, 255], "gray51"=>[130, 130, 130], "rosybrown4"=>[139, 105, 105], "sienna4"=>[139, 71, 38], "paleturquoise3"=>[150, 205, 205], "green2"=>[0, 238, 0], "pink3"=>[205, 145, 158], "wheat4"=>[139, 126, 102], "lightskyblue2"=>[164, 211, 238], "yellow1"=>[255, 255, 0], "sgilightgray"=>[170, 170, 170], "paleturquoise4"=>[102, 139, 139], "green3"=>[0, 205, 0], "pink4"=>[139, 99, 108], "lightskyblue3"=>[141, 182, 205], "yellow2"=>[238, 238, 0], "gray52"=>[133, 133, 133], "mediumorchid"=>[186, 85, 211], "green4"=>[0, 139, 0], "lightskyblue4"=>[96, 123, 139], "ivory1"=>[255, 255, 240], "yellow3"=>[205, 205, 0], "plum1"=>[255, 187, 255], "lightblue1"=>[191, 239, 255], "seashell1"=>[255, 245, 238], "gray53"=>[135, 135, 135], "ivory2"=>[238, 238, 224], "violetred1"=>[255, 62, 150], "plum2"=>[238, 174, 238], "lightblue2"=>[178, 223, 238], "seashell2"=>[238, 229, 222], "gray54"=>[138, 138, 138], "white"=>[255, 255, 255], "mediumslateblue"=>[123, 104, 238], "ivory3"=>[205, 205, 193], "yellow4"=>[139, 139, 0], "lightcyan1"=>[224, 255, 255], "plum3"=>[205, 150, 205], "lightslateblue"=>[132, 112, 255], "lightblue3"=>[154, 192, 205], "seashell3"=>[205, 197, 191], "gray55"=>[140, 140, 140], "ivory4"=>[139, 139, 131], "lightcyan2"=>[209, 238, 238], "violetred2"=>[238, 58, 140], "peru"=>[205, 133, 63], "lightblue4"=>[104, 131, 139], "seashell4"=>[139, 134, 130], "gray56"=>[143, 143, 143], "flesh"=>[255, 125, 64], "lightcyan3"=>[180, 205, 205], "violetred3"=>[205, 50, 120], "plum4"=>[139, 102, 139], "lightblue"=>[173, 216, 230], "gray57"=>[145, 145, 145], "coral"=>[255, 127, 80], "cadmiumorange"=>[255, 97, 3], "khaki1"=>[255, 246, 143], "lightcyan4"=>[122, 139, 139], "violetred4"=>[139, 34, 82], "lightslategray"=>[119, 136, 153], "gray58"=>[148, 148, 148], "sgigray32"=>[81, 81, 81], "khaki2"=>[238, 230, 133], "dodgerblue1"=>[30, 144, 255], "lavenderblush1"=>[255, 240, 245], "gray59"=>[150, 150, 150], "khaki3"=>[205, 198, 115], "dodgerblue2"=>[28, 134, 238], "lavenderblush2"=>[238, 224, 229], "banana"=>[227, 207, 87], "khaki4"=>[139, 134, 78], "dodgerblue3"=>[24, 116, 205], "lavenderblush3"=>[205, 193, 197], "cyan2"=>[0, 238, 238], "cadetblue"=>[95, 158, 160], "cadmiumyellow"=>[255, 153, 18], "mediumpurple"=>[147, 112, 219], "dodgerblue4"=>[16, 78, 139], "lavenderblush4"=>[139, 131, 134], "cyan3"=>[0, 205, 205], "silver"=>[192, 192, 192], "sgigray36"=>[91, 91, 91], "cyan4"=>[0, 139, 139], "deeppink1"=>[255, 20, 147], "burntumber"=>[138, 51, 36], "deeppink2"=>[238, 18, 137], "aliceblue"=>[240, 248, 255], "rosybrown"=>[188, 143, 143], "aquamarine"=>[127, 255, 212], "plum"=>[221, 160, 221], "deeppink3"=>[205, 16, 118], "darkgoldenrod"=>[184, 134, 11], "darkslategray1"=>[151, 255, 255], "chartreuse1"=>[127, 255, 0], "deeppink4"=>[139, 10, 80], "antiquewhite"=>[250, 235, 215], "darkslategray2"=>[141, 238, 238], "chartreuse2"=>[118, 238, 0], "darkslategray3"=>[121, 205, 205], "raspberry"=>[135, 38, 87], "chartreuse3"=>[102, 205, 0], "darkslategray4"=>[82, 139, 139], "chartreuse4"=>[69, 139, 0], "lightsalmon"=>[255, 160, 122], "warmgrey"=>[128, 128, 105], "cobalt"=>[61, 89, 171], "pink"=>[255, 192, 203], "orangered"=>[255, 69, 0], "indianred"=>[205, 92, 92], "palegreen1"=>[154, 255, 154], "red1"=>[255, 0, 0], "darkorange1"=>[255, 127, 0], "palegreen2"=>[144, 238, 144], "darkorchid"=>[153, 50, 204], "red2"=>[238, 0, 0], "darkorange2"=>[238, 118, 0], "palegreen3"=>[124, 205, 124], "goldenrod"=>[218, 165, 32], "red3"=>[205, 0, 0], "tomato1"=>[255, 99, 71], "yellowgreen"=>[154, 205, 50], "darkorange3"=>[205, 102, 0], "palegreen4"=>[84, 139, 84], "gray30"=>[77, 77, 77], "red4"=>[139, 0, 0], "tomato2"=>[238, 92, 66], "darkorange4"=>[139, 69, 0], "antiquewhite1"=>[255, 239, 219], "gray31"=>[79, 79, 79], "brown"=>[165, 42, 42], "tomato3"=>[205, 79, 57], "thistle1"=>[255, 225, 255], "antiquewhite2"=>[238, 223, 204], "darkolivegreen1"=>[202, 255, 112], "gray32"=>[82, 82, 82], "gray80"=>[204, 204, 204], "tomato4"=>[139, 54, 38], "sepia"=>[94, 38, 18], "thistle2"=>[238, 210, 238], "deepskyblue"=>[0, 191, 255], "antiquewhite3"=>[205, 192, 176], "darkolivegreen2"=>[188, 238, 104], "cornflowerblue"=>[100, 149, 237], "gray33"=>[84, 84, 84], "gray81"=>[207, 207, 207], "thistle3"=>[205, 181, 205], "peachpuff"=>[255, 218, 185], "steelblue1"=>[99, 184, 255], "antiquewhite4"=>[139, 131, 120], "darkolivegreen3"=>[162, 205, 90], "gray34"=>[87, 87, 87], "gray82"=>[209, 209, 209], "maroon"=>[128, 0, 0], "thistle4"=>[139, 123, 139], "springgreen1"=>[0, 238, 118], "steelblue2"=>[92, 172, 238], "darkolivegreen4"=>[110, 139, 61], "gray35"=>[89, 89, 89], "gray83"=>[212, 212, 212], "springgreen2"=>[0, 205, 102], "steelblue3"=>[79, 148, 205], "steelblue"=>[70, 130, 180], "gray36"=>[92, 92, 92], "gray84"=>[214, 214, 214], "springgreen3"=>[0, 139, 69], "steelblue4"=>[54, 100, 139], "gray37"=>[94, 94, 94], "blue2"=>[0, 0, 238], "darkturquoise"=>[0, 206, 209], "gray38"=>[97, 97, 97], "gray85"=>[217, 217, 217], "sgigray12"=>[30, 30, 30], "blanchedalmond"=>[255, 235, 205], "violet"=>[238, 130, 238], "blue3"=>[0, 0, 205], "gray39"=>[99, 99, 99], "gray86"=>[219, 219, 219], "burlywood1"=>[255, 211, 155], "navy"=>[0, 0, 128], "darkviolet"=>[148, 0, 211], "blue4"=>[0, 0, 139], "gray87"=>[222, 222, 222], "burlywood2"=>[238, 197, 145], "floralwhite"=>[255, 250, 240], "turquoise"=>[64, 224, 208], "darkorange"=>[255, 140, 0], "gray88"=>[224, 224, 224], "burlywood3"=>[205, 170, 125], "gray89"=>[227, 227, 227], "dimgray"=>[105, 105, 105], "sgigray16"=>[40, 40, 40], "snow"=>[255, 250, 250], "burlywood4"=>[139, 115, 85], "mediumvioletred"=>[199, 21, 133], "seagreen"=>[46, 139, 87], "purple1"=>[155, 48, 255], "lightsteelblue1"=>[202, 225, 255], "purple2"=>[145, 44, 238], "orangered1"=>[255, 69, 0], "lightsteelblue2"=>[188, 210, 238], "purple3"=>[125, 38, 205], "cyan/aqua"=>[0, 255, 255], "orangered2"=>[238, 64, 0], "lightsteelblue3"=>[162, 181, 205], "purple4"=>[85, 26, 139], "orangered3"=>[205, 55, 0], "fuchsia"=>[255, 0, 255], "lightsteelblue4"=>[110, 123, 139], "mistyrose1"=>[255, 228, 225], "orangered4"=>[139, 37, 0], "darkorchid1"=>[191, 62, 255], "lemonchiffon1"=>[255, 250, 205], "mistyrose2"=>[238, 213, 210], "darkorchid2"=>[178, 58, 238], "lemonchiffon2"=>[238, 233, 191], "mistyrose3"=>[205, 183, 181], "royalblue"=>[65, 105, 225], "darkorchid3"=>[154, 50, 205], "lemonchiffon3"=>[205, 201, 165], "sgiolivedrab"=>[142, 142, 56], "mistyrose4"=>[139, 125, 123], "darkorchid4"=>[104, 34, 139], "peachpuff1"=>[255, 218, 185], "mediumseagreen"=>[60, 179, 113], "lemonchiffon4"=>[139, 137, 112], "lightgreen"=>[144, 238, 144], "peachpuff2"=>[238, 203, 173], "peachpuff3"=>[205, 175, 149], "peachpuff4"=>[139, 119, 101], "lightsteelblue"=>[176, 196, 222], "skyblue"=>[135, 206, 235], "gray10"=>[26, 26, 26], "magenta"=>[255, 0, 255], "turquoise1"=>[0, 245, 255], "gray11"=>[28, 28, 28], "turquoise2"=>[0, 229, 238], "chartreuse"=>[127, 255, 0], "navajowhite"=>[255, 222, 173], "ghostwhite"=>[248, 248, 255], "gray12"=>[31, 31, 31], "gray60"=>[153, 153, 153], "melon"=>[227, 168, 105], "turquoise3"=>[0, 197, 205], "gray13"=>[33, 33, 33], "gray61"=>[156, 156, 156], "coral1"=>[255, 114, 86], "turquoise4"=>[0, 134, 139], "gray14"=>[36, 36, 36], "gray62"=>[158, 158, 158], "sgiteal"=>[56, 142, 142], "coral2"=>[238, 106, 80], "honeydew"=>[240, 255, 240], "gray1"=>[3, 3, 3], "gray15"=>[38, 38, 38], "firebrick"=>[178, 34, 34], "coral3"=>[205, 91, 69], "lightpink"=>[255, 182, 193], "gray2"=>[5, 5, 5], "gray16"=>[41, 41, 41], "gray63"=>[161, 161, 161], "coral4"=>[139, 62, 47], "mediumblue"=>[0, 0, 205], "gray3"=>[8, 8, 8], "gray17"=>[43, 43, 43], "gray64"=>[163, 163, 163], "blueviolet"=>[138, 43, 226], "darkgreen"=>[0, 100, 0], "gray4"=>[10, 10, 10], "gray18"=>[46, 46, 46], "gray65"=>[166, 166, 166], "deepskyblue1"=>[0, 191, 255], "gray5"=>[13, 13, 13], "gray19"=>[48, 48, 48], "gray66"=>[168, 168, 168], "mediumspringgreen"=>[0, 250, 154], "darkkhaki"=>[189, 183, 107], "lightgoldenrodyellow"=>[250, 250, 210], "orange"=>[255, 128, 0], "deepskyblue2"=>[0, 178, 238], "gray6"=>[15, 15, 15], "gray67"=>[171, 171, 171], "mediumturquoise"=>[72, 209, 204], "deepskyblue3"=>[0, 154, 205], "turquoiseblue"=>[0, 199, 140], "gray7"=>[18, 18, 18], "gray68"=>[173, 173, 173], "orchid"=>[218, 112, 214], "burlywood"=>[222, 184, 135], "thistle"=>[216, 191, 216], "mediumpurple1"=>[171, 130, 255], "deepskyblue4"=>[0, 104, 139], "gray8"=>[20, 20, 20], "gray69"=>[176, 176, 176], "palevioletred"=>[219, 112, 147], "sandybrown"=>[244, 164, 96], "mediumpurple2"=>[159, 121, 238], "yellow"=>[255, 255, 0], "darkgoldenrod1"=>[255, 185, 15], "gray9"=>[23, 23, 23], "darkseagreen"=>[143, 188, 143], "mediumpurple3"=>[137, 104, 205], "darkgoldenrod2"=>[238, 173, 14], "sgigray92"=>[234, 234, 234], "mediumpurple4"=>[93, 71, 139], "darkgoldenrod3"=>[205, 149, 12], "mediumaquamarine"=>[102, 205, 170], "firebrick1"=>[255, 48, 48], "tomato"=>[255, 99, 71], "darkgoldenrod4"=>[139, 101, 8], "gray"=>[128, 128, 128], "firebrick2"=>[238, 44, 44], "palegoldenrod"=>[238, 232, 170], "moccasin"=>[255, 228, 181], "firebrick3"=>[205, 38, 38], "snow1"=>[255, 250, 250], "sgigray96"=>[244, 244, 244], "firebrick4"=>[139, 26, 26], "snow2"=>[238, 233, 233], "papayawhip"=>[255, 239, 213], "snow3"=>[205, 201, 201], "eggshell"=>[252, 230, 201], "snow4"=>[139, 137, 137], "tan"=>[210, 180, 140], "peacock"=>[51, 161, 201], "gold1"=>[255, 215, 0], "slateblue1"=>[131, 111, 255], "springgreen"=>[0, 255, 127], "ivoryblack"=>[41, 36, 33], "gold2"=>[238, 201, 0], "honeydew1"=>[240, 255, 240], "greenyellow"=>[173, 255, 47], "slateblue2"=>[122, 103, 238], "gainsboro"=>[220, 220, 220], "azure"=>[240, 255, 255], "gold3"=>[205, 173, 0], "oldlace"=>[253, 245, 230], "honeydew2"=>[224, 238, 224], "gold4"=>[139, 117, 0], "honeydew3"=>[193, 205, 193], "slateblue3"=>[105, 89, 205], "green"=>[0, 128, 0], "skyblue1"=>[135, 206, 255], "honeydew4"=>[131, 139, 131], "cadetblue1"=>[152, 245, 255], "slateblue4"=>[71, 60, 139], "sienna"=>[160, 82, 45], "darksalmon"=>[233, 150, 122], "skyblue2"=>[126, 192, 238], "cadetblue2"=>[142, 229, 238], "sapgreen"=>[48, 128, 20], "gray40"=>[102, 102, 102], "skyblue3"=>[108, 166, 205], "midnightblue"=>[25, 25, 112], "blue"=>[0, 0, 255], "cadetblue3"=>[122, 197, 205], "sgichartreuse"=>[113, 198, 113], "skyblue4"=>[74, 112, 139], "cadetblue4"=>[83, 134, 139], "lime"=>[0, 255, 0], "gray90"=>[229, 229, 229], "gray42"=>[105, 105, 105], "mistyrose"=>[255, 228, 225], "burntsienna"=>[138, 54, 15], "ivory"=>[255, 255, 240], "gray43"=>[110, 110, 110], "gray91"=>[232, 232, 232], "darkcyan"=>[0, 139, 139], "gray44"=>[112, 112, 112], "gray92"=>[235, 235, 235], "lawngreen"=>[124, 252, 0], "gray45"=>[115, 115, 115], "gray93"=>[237, 237, 237], "palevioletred1"=>[255, 130, 171], "saddlebrown"=>[139, 69, 19], "lavenderblush"=>[255, 240, 245], "gray46"=>[117, 117, 117], "gray94"=>[240, 240, 240], "palevioletred2"=>[238, 121, 159], "khaki"=>[240, 230, 140], "paleturquoise"=>[174, 238, 238], "gray47"=>[120, 120, 120], "gray95"=>[242, 242, 242], "carrot"=>[237, 145, 33], "palevioletred3"=>[205, 104, 137], "lightsalmon1"=>[255, 160, 122], "olivedrab"=>[107, 142, 35], "gray48"=>[122, 122, 122], "palevioletred4"=>[139, 71, 93], "wheat"=>[245, 222, 179], "lightsalmon2"=>[238, 149, 114], "powderblue"=>[176, 224, 230], "gray49"=>[125, 125, 125], "gray96"=>[245, 245, 245], "red"=>[255, 0, 0], "brown1"=>[255, 64, 64], "lightsalmon3"=>[205, 129, 98], "manganeseblue"=>[3, 168, 158], "gray97"=>[247, 247, 247], "brown2"=>[238, 59, 59], "lightsalmon4"=>[139, 87, 66], "olive"=>[128, 128, 0], "gray98"=>[250, 250, 250], "lightgrey"=>[211, 211, 211], "sgigray72"=>[183, 183, 183], "brown3"=>[205, 51, 51], "gold"=>[255, 215, 0], "mintcream"=>[245, 255, 250], "gray99"=>[252, 252, 252], "brown4"=>[139, 35, 35], "sgidarkgray"=>[85, 85, 85], "magenta2"=>[238, 0, 238], "sgigray76"=>[193, 193, 193], "magenta3"=>[205, 0, 205], "dodgerblue"=>[30, 144, 255], "magenta4"=>[139, 0, 139], "cobaltgreen"=>[61, 145, 64], "sgisalmon"=>[198, 113, 113], "emeraldgreen"=>[0, 201, 87], "sgibrightgray"=>[197, 193, 170], "darkred"=>[139, 0, 0], "rawsienna"=>[199, 97, 20], "seagreen1"=>[84, 255, 159], "lightyellow"=>[255, 255, 224], "olivedrab1"=>[192, 255, 62], "orchid1"=>[255, 131, 250], "seashell"=>[255, 245, 238], "seagreen2"=>[78, 238, 148], "bisque1"=>[255, 228, 196], "olivedrab2"=>[179, 238, 58], "orchid2"=>[238, 122, 233], "azure1"=>[240, 255, 255], "seagreen3"=>[67, 205, 128], "bisque2"=>[238, 213, 183], "whitesmoke"=>[245, 245, 245], "olivedrab3"=>[154, 205, 50], "orchid3"=>[205, 105, 201], "azure2"=>[224, 238, 238], "lightpink1"=>[255, 174, 185], "lightgoldenrod1"=>[255, 236, 139], "bisque3"=>[205, 183, 158], "olivedrab4"=>[105, 139, 34], "orchid4"=>[139, 71, 137], "purple"=>[128, 0, 128], "lemonchiffon"=>[255, 250, 205], "azure3"=>[193, 205, 205], "lightpink2"=>[238, 162, 173], "seagreen4"=>[46, 139, 87], "forestgreen"=>[34, 139, 34], "lightgoldenrod2"=>[238, 220, 130], "bisque4"=>[139, 125, 107], "sgislateblue"=>[113, 113, 198], "azure4"=>[131, 139, 139], "lightpink3"=>[205, 140, 149], "lightgoldenrod3"=>[205, 190, 112], "gray20"=>[51, 51, 51], "sgibeet"=>[142, 56, 142], "limegreen"=>[50, 205, 50], "orange1"=>[255, 165, 0], "lightpink4"=>[139, 95, 101], "lightgoldenrod4"=>[139, 129, 76], "gray21"=>[54, 54, 54], "orange2"=>[238, 154, 0], "gray22"=>[56, 56, 56], "gray70"=>[179, 179, 179], "orange3"=>[205, 133, 0], "gray23"=>[59, 59, 59], "gray71"=>[181, 181, 181], "orange4"=>[139, 90, 0], "lightseagreen"=>[32, 178, 170], "bisque"=>[255, 228, 196], "gray24"=>[61, 61, 61], "gray72"=>[184, 184, 184], "darkslateblue"=>[72, 61, 139], "gray25"=>[64, 64, 64], "gray73"=>[186, 186, 186], "lightcoral"=>[240, 128, 128], "tan1"=>[255, 165, 79], "gray26"=>[66, 66, 66], "black"=>[0, 0, 0], "tan2"=>[238, 154, 73], "gray27"=>[69, 69, 69], "gray74"=>[189, 189, 189], "teal"=>[0, 128, 128], "tan3"=>[205, 133, 63], "aquamarine1"=>[127, 255, 212], "gray28"=>[71, 71, 71], "gray75"=>[191, 191, 191], "tan4"=>[139, 90, 43], "darkmagenta"=>[139, 0, 139], "aquamarine2"=>[118, 238, 198], "gray29"=>[74, 74, 74], "gray76"=>[194, 194, 194], "indianred1"=>[255, 106, 106], "chocolate"=>[210, 105, 30], "aquamarine3"=>[102, 205, 170], "gray77"=>[196, 196, 196], "indianred2"=>[238, 99, 99], "slategray1"=>[198, 226, 255], "chocolate1"=>[255, 127, 36], "aquamarine4"=>[69, 139, 116], "gray78"=>[199, 199, 199], "sgigray52"=>[132, 132, 132], "indianred3"=>[205, 85, 85], "salmon"=>[250, 128, 114], "darkblue"=>[0, 0, 139], "lightyellow1"=>[255, 255, 224], "slategray2"=>[185, 211, 238], "chocolate2"=>[238, 118, 33], "gray79"=>[201, 201, 201], "indianred4"=>[139, 58, 58], "hotpink"=>[255, 105, 180], "brick"=>[156, 102, 31], "lightyellow2"=>[238, 238, 209], "slategray3"=>[159, 182, 205], "chocolate3"=>[205, 102, 29], "darkgray"=>[169, 169, 169], "lightyellow3"=>[205, 205, 180], "darkolivegreen"=>[85, 107, 47], "slategray4"=>[108, 123, 139], "chocolate4"=>[139, 69, 19], "salmon1"=>[255, 140, 105], "deeppink"=>[255, 20, 147], "indigo"=>[75, 0, 130], "darkslategray"=>[47, 79, 79], "navajowhite1"=>[255, 222, 173], "lightyellow4"=>[139, 139, 122], "royalblue1"=>[72, 118, 255], "slateblue"=>[106, 90, 205], "maroon1"=>[255, 52, 179], "sgigray56"=>[142, 142, 142], "salmon2"=>[238, 130, 98], "navajowhite2"=>[238, 207, 161], "coldgrey"=>[128, 138, 135], "royalblue2"=>[67, 110, 238], "cornsilk1"=>[255, 248, 220], "maroon2"=>[238, 48, 167], "salmon3"=>[205, 112, 84], "navajowhite3"=>[205, 179, 139], "royalblue3"=>[58, 95, 205], "cornsilk2"=>[238, 232, 205], "maroon3"=>[205, 41, 144], "slategray"=>[112, 128, 144], "salmon4"=>[139, 76, 57], "navajowhite4"=>[139, 121, 94], "royalblue4"=>[39, 64, 139], "darkseagreen1"=>[193, 255, 193], "cornsilk3"=>[205, 200, 177], "maroon4"=>[139, 28, 98], "hotpink1"=>[255, 110, 180], "mediumorchid1"=>[224, 102, 255], "darkseagreen2"=>[180, 238, 180], "cornsilk4"=>[139, 136, 120], "mint"=>[189, 252, 201], "hotpink2"=>[238, 106, 167], "sgilightblue"=>[125, 158, 192], "goldenrod1"=>[255, 193, 37], "mediumorchid2"=>[209, 95, 238], "darkseagreen3"=>[155, 205, 155], "hotpink3"=>[205, 96, 144], "goldenrod2"=>[238, 180, 34], "mediumorchid3"=>[180, 82, 205], "darkseagreen4"=>[105, 139, 105], "hotpink4"=>[139, 58, 98], "crimson"=>[220, 20, 60], "lightskyblue"=>[135, 206, 250], "goldenrod3"=>[205, 155, 29], "mediumorchid4"=>[122, 55, 139], "cornsilk"=>[255, 248, 220]}
        pdf=FPDF.new('P', "cm")
        pdf.SetTitle(pdfFile)
        pdf.SetCreator("newickDraw")
        pdf.SetAuthor(ENV["USER"]) if (!ENV["USER"].nil?)
        pdf.AddPage
        yUnit = nil
        lineWidth = nil
        fontSize = nil
        bootScale = 0.6
        if (taxa.size < 30)
          fontSize = 10
          yUnit = 0.5
          lineWidth = 0.02
        elsif (taxa.size < 60)
          fontSize = 8
          yUnit = 0.25
          lineWidth = 0.01
        elsif (taxa.size < 150)
          fontSize = 6
          yUnit = 0.197
          lineWidth = 0.01
        elsif (taxa.size < 300)
          fontSize = 2
          yUnit = 0.09
          lineWidth = 0.005
        elsif (taxa.size < 400)
          fontSize = 2
          yUnit = 0.055
          lineWidth = 0.002
        elsif (taxa.size < 800)
          fontSize = 1
          yUnit = 0.030
          lineWidth = 0.0015
        else
          fontSize = 0.5
          yUnit = 0.020
          lineWidth = 0.0010
        end
        bootScale = 0.5 * fontSize
        pdf.SetFont('Times','B', fontSize)
        calcPos(yUnit) # calculate node pos before drawing
        max = 0
        @root.leaves.each do |leaf|
          d = leaf.distToAncestor(@root)
          max = d if (max < d)
        end
        xScale = 10.0/max
        xOffSet = 0.25
        pdf.SetLineWidth(lineWidth)
        pdf.SetTextColor(0, 0, 0)
        pdf.Line(0, @root.y, xOffSet, @root.y)
        pdf.Line(xOffSet, @root.yMin, xOffSet, @root.yMax)
        @root.descendants.each do |child|
          if (!child.leaf?)
            if (child.name.to_f >= goodBoot && boot == "width") # good bootstrap
              pdf.SetLineWidth(lineWidth * 5)
            else
              pdf.SetLineWidth(lineWidth)
            end
            bootX = xOffSet + child.x*xScale 
            bootY = ((child.yMin + child.yMax) / 2.0) 
            pdf.SetXY(bootX, bootY) 
            pdf.SetFont('Times','B', bootScale)
            pdf.Write(0, child.name.to_s)
            pdf.SetFont('Times','B', fontSize)
            pdf.Line(xOffSet + child.x*xScale, child.yMin, 
            xOffSet + child.x*xScale, child.yMax)
          else
            if (child.parent.name.to_f > goodBoot && boot == "width") # good bootstrap
              pdf.SetLineWidth(lineWidth * 5)
            else
              pdf.SetLineWidth(lineWidth)
            end
            pdf.SetXY(xOffSet + child.x*xScale, child.y)
            efields = child.name.split("__")
            entry, species = efields.first, efields.last
            if (entry =~/\{([^\}]*)\}/)
              species = $1
            end
            species = entry if species.nil? && !rawNames
            species = child.name if rawNames
            hl = false
            highlights.keys.each do |highlight|
              hl = highlights[highlight] if (entry.index(highlight))
            end
            if (pdfFile.index(entry)) # name of query taxon
              pdf.SetTextColor(colors["red"][0], colors["red"][1], colors["red"][2]) 
              pdf.Write(0, entry) 
              pdf.SetTextColor(colors["black"][0], colors["black"][1], colors["black"][2]) 
            elsif (linker && link = send(linker, entry)) 
              pdf.SetTextColor(colors["black"][0], colors["black"][1], colors["black"][2]) if hl 
              pdf.Write(0, species, link)
              pdf.SetTextColor(colors["black"][0], colors["black"][1], colors["black"][2]) if hl 
            elsif (!species.nil?)
              pdf.SetTextColor(colors[hl][0], colors[hl][1], colors[hl][2]) if hl 
              pdf.Write(0, species) 
              pdf.SetTextColor(colors["black"][0], colors["black"][1], colors["black"][2]) if hl 
            else
              pdf.SetTextColor(colors[hl][0], colors[hl][1], colors[hl][2]) if hl 
              pdf.Write(0, entry) 
              pdf.SetTextColor(colors["black"][0], colors["black"][1], colors["black"][2]) if hl  
            end
          end
          pdf.Line(xOffSet + child.parent.x*xScale, child.y, 
          xOffSet + child.x*xScale, child.y)
        end
        if (labelName)
          pdf.SetFont('Times','B', 24)
          pdf.SetXY(0, pdf.GetY + 1)
          pdf.Write(0, labelName)
        end
        if (brackets)
          brackets.each do |bracket|
            x, y1, y2, label, r, p = bracket
            next if label.nil?
            pdf.SetLineWidth(lineWidth * 5)
            pdf.SetFont('Times','B', fontSize*1.5)
            pdf.Line(x, y1, x, y2)
            pdf.Line(x, y1, x - 0.3, y1)
            pdf.Line(x, y2, x - 0.3, y2)
            pdf.SetXY(x, (y1+y2)/2)
            pdf.Write(0, label)
            if (r == "r")
              pdf.SetTextColor(255, 0, 0) 
              pdf.SetXY(x + 1.8, -0.65+(y1+y2)/2)
              pdf.SetFont('Times','B', fontSize*10)
              pdf.Write(0, " .")
              pdf.SetTextColor(0, 0, 0)
            end
            if (p == "p" || r == "p")
              pdf.SetTextColor(255, 0, 255) 
              pdf.SetXY(x + 2.3, -0.65+(y1+y2)/2)
              pdf.SetFont('Times','B', fontSize*10)
              pdf.Write(0, " .")
              pdf.SetTextColor(0, 0, 0)
            end
          end
        end
        pdf.SetLineWidth(lineWidth * 5)
        pdf.Line(1, pdf.GetY + 1, 1 + 0.1*xScale, pdf.GetY + 1) 
        pdf.SetFont('Times','B', fontSize)
        pdf.SetXY(1 + 0.1*xScale, pdf.GetY + 1)
        pdf.Write(0, "0.1")
        if (pdfFile =~/^--/)
          return pdf.Output
        else
          pdf.Output(pdfFile)
        end
      end
    end



