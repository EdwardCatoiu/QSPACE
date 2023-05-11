import numpy as np

def bfs_algo(parent, children_list,answers =[], incomplete_matches = [], over_matches = [], enzyme =''):
    
    def buildLevel(queue, children_list):
        """
        #title: dynamic BFS with no redundancy
        #author: Maxwell Lu & Edward Catoiu
        #objective: take multiple children and find all combinations of children to obtain the parent node of key:value pairs


        """
        global total_alive_leaves
        levelHash = {}
        level = []
        node_count = 0
        leaf_count = 0
        for node in queue:
#             print '\n NODE:' ,node_count , node
            node_count +=1 
            #create level by adding children to each queued item
            for child in children_list:
                #checks go here
#                 print ' \t CHILD:', child
                leaf = add(node, child)
                leaf_count += 1
                leafHash = buildHash(leaf, children_list)
                #check if existis else add
                #print leafHash , 'leafhash'
                
#                 print ' \t LEAF: ', leaf_count,  leaf
               # print levelHash, 'levelhash'
                if leafHash not in levelHash:
                    #add leafHash & exec code
                    levelHash[leafHash] = 1
#                     print leaf, 'leaf'
                    if alive(leaf,parent, leafHash):
                        level.append(leaf)
                        #these are the combinations that have not quite been killed
#                         print 'alive LEAF -->' ,leafHash
#                         total_alive_leaves.append(leafHash)
#                     else: 
#                         print 'leaf died here -->' , leafHash
                   

        # print level;
        return level;

    def buildHash(node,children_list):
#         print 'node:' , node
        children = list(np.zeros(len(children_list),dtype = int))

        #children = [0, 0, 0, 0, 0, 0]
        for child in node['children']:
            children[child['child']] += 1
            # print ''.join(map(lambda x: str(x),children))
        return ''.join(map(lambda x: str(x),children))

    #eddie's main function
    #designed to be used for generalized "parent" dictionary 
    def main(parent,children_list):
        temp = {}
        root = {}
        for key in parent.keys():
            temp.update({key : 0})
        root.update({'sum' : temp})
        root.update({'children' : []})

        queue = [root]
        while len(queue) > 0:
            #process queue
            queue = buildLevel(queue,children_list)

    #eddie's sum function, generalized to any key:value pair
    def add(a, b):
        temp = {}
        for key in a['sum'].keys():
            sum_value = a['sum'][key] + b[key]
            temp.update({key : sum_value})
            
        for key in b:
            if key not in a['sum'].keys() and key != 'child':
                temp.update({key : b[key]})

        return {
            'sum' : temp,
            'children' : a['children'] + [b]
        }

    #eddie's alive node
    def alive(node, parent,leafHash):
        #global answers
        # compare to Parent
        #print node['sum']
        if node['sum'] == parent:
            answers.append(node['children'])
            return False

        for key in node['sum']:
            if key not in parent.keys():
                over_matches.append(leafHash)
                return False
            
            elif node['sum'][key] > parent[key]:
                over_matches.append(leafHash)
                return False
        incomplete_matches.append(leafHash)
        return True
    
    main(parent,children_list)
    
    answers_string =[]
    for ans in answers:
        children = list(np.zeros(len(children_list),dtype = int))
        for child in ans:
            children[child['child']] += 1
        answers_string.append('_'.join(map(lambda x: str(x),children)))
    #IF THERE ARE FULL MATCHES, just return them!!!!!
    if answers_string !=[]:
        return answers,answers_string, [], []
    
    def filter_over_matches(over_matches,e=enzyme):
        """
        finds over matches which only use 1 pdb. these are biologically relevant
        """
        biological_over_matches = []
        printing = False
        if len(over_matches) >700:
            print ('\n' , e)
            printing = True
            print ('Before Filtering Overmatches: %i' %len(over_matches))
        
        for o in over_matches:
            #if np.sum([int(s) for s in o.split('_')]) > 1:
                #continue
                
            if o.count('1') != 1 or o.count('0') != len(o) - 1:
                continue
            biological_over_matches.append(o)
        if printing:
            print ('After Filtering Overmatches: %i' %len(biological_over_matches))
        return biological_over_matches

    over_string = []
    if over_matches != []:
        over_matches = filter_over_matches(over_matches , e = enzyme)
        
         #o -> to string
        for o_match in over_matches:
            new_o_match= ''
            for i, val in enumerate(o_match):
                new_o_match = new_o_match + val + '_'

            new_o_match = new_o_match.rstrip('_')
            #print new_i_match
            over_string.append(new_o_match)
   
        
    incomplete_string = []   
    if incomplete_matches != []:
        incomplete_matches = filter_incomplete_matches(incomplete_matches,e=enzyme)

        #i -> to string
        
        for i_match in incomplete_matches:
            new_i_match= ''
            for i, val in enumerate(i_match):
                new_i_match = new_i_match + val + '_'

            new_i_match = new_i_match.rstrip('_')
            #print new_i_match
            incomplete_string.append(new_i_match)
    
   
        
    return answers,answers_string, incomplete_string, over_string