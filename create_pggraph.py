


	def parseFeatures(self, feature_file):
		print "parsing features and constructing kmer graph"	
		ihandle=open(feature_file, 'r')
		kmer_q=deque()
		prev_seq=""
		cur_seq=""
		prev_kmer=None
		#loop through figfams to create kmers
		for line in ihandle:
			if line.startswith('#'):
				continue

			feature=geneInfo(line=line)
            if org_id not in self.replicon_map:
				self.replicon_map[org_id]=set()
			else:
				self.replicon_map[org_id].add(replicon_id)

			if ignore and fam_id in self.ignore_fams:
				continue
			num_fam+=1
			cur_seq=replicon_id
			self.replicon_ids.add(cur_seq)
			if(prev_seq != cur_seq and num_fam>1):
				kmer_q=deque()#clear kmer stack because switching replicons
				prev_kmer=None

			kmer_q.append(fig_info)#append the figfam ID
			if(len(kmer_q)>self.ksize):
				kmer_q.popleft()
				kmer=self.addKmer(prev_kmer, kmer_q)
			elif(len(kmer_q)== self.ksize):
				kmer=self.addKmer(prev_kmer, kmer_q)#right now only passing in the last figfams information
			else:#kmer size is less than ksize
				kmer=None
			prev_seq=cur_seq # record which replicon we are on
			prev_kmer=kmer
		inHandle.close()

    def hashKmer(self, feature_list):
        reverese, palindrome, feature_list = flipKmer(feature_list)
        kmer_hash=[]
        for feature in feature_list:
            if not feature.group_name in self.groups_seen:
                feature.group_id=len(self.group_index)
                self.groups_seen[feature.group_name]=feature.group_id
                self.group_index.append(feature.group_name)
            else:
                feature.group_id=self.groups_seen[feature.group_name]
            kmer_hash.append(feature.group_id)
            feature.feature_id= len(self.feature_index)
            self.feature_index.append(feature)
        return kmer_hash
            


    def flipKmer(self, feature_list)
        i=0
        k_size=len(id_list)
        palindrome=0
        reverse=0
        while i<(len(k_size)/2):
            if id_list[i].group_name < id_list[k_size-(i+1)].group_name:
                return (reverse, palidrome, feature_list)
            else if feature_list[i].group_name > feature_list[k_size-(i+1)].group_name:
                reverse=1
                return (reverse, palindrome, feature_list.reverse())
            i+=1
        palindrome =1
        return (reverse, palindrome, feature_list)


def RF_to_PG(rfgraph, non_dups):
    for cur_node in non_dups:
        DFSExpand(cur_node, rfgraph)



def DFSExpand(cur_node, rfgraph, thread_bundle):
    if cur_node.anchor:

    #Incoming
    #for existing 'threads' gather edges, assign to pg-nodes, add to queue
    if consistent and len(prev_assigned.forward) >0:
        adjust_by=range(1,k+1)
    for feat in prev_assigned.forward:
        for i in adjust_by:
            if feat+i in cur_node.family_members[i-1]:
                #make assigned
        
    visit_queue=set()
    if not cur_node.duplicate:
        storage.ExpandNode(cur_node)
    for fam in cur_node.family_members:
        for feat in fam:
            visit_queue.add(feat)
    #Processing

    #Outgoing
    #Receiving
    #Returning
    
