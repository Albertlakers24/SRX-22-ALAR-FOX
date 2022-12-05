import numpy as np
#import matplotlib.pyplot as plt

dict = {
    "subsystem": ["jet fuel", "SAF", "Liquid Hydrogen", "Electric"],
    "criteria" : ["Unit/prod. cost", "Energy volume density", "Energy mass density", "Maturity", "Integration feasibility","Emissions"],

}

weights = np.arange(1,10,1)

rating = np.transpose(np.array(([6,5,9,3],
                 [9,9,5,3],
                 [6,6,8,2],
                 [9,8,5,5],
                 [9,8,3,4],
                 [2,5,8,8])))
original_rating = np.array((7,8,8,6,6,9))
#print(rating[4][0])
# criteria = dict["criteria"][0]
# subsystem = dict["subsystem"][0]

# original_value_list =[[42,35,63, 21],
#                       [72,72,40,24],
#                       [48,48,64,16],
#                       [54,48,30,30],
#                       [54,48,18,24],
#                       [18,45,72,72]]
original_value_list = np.transpose(original_rating*rating)

original_sum = np.sum(original_value_list,axis=0)

#print(original_value_list)

# for i in range(len(dict["subsystem"])):
#     for j in range(len(dict["criteria"])):
#         subsystem = dict["subsystem"][i]
#         criteria = dict["criteria"][j]
#         value = rating[j][i] * original_rating[j]
#         original_value_list.append(value)
# print(original_value_list)

value_list = [[[None]*9]*4]*6

# for j in range(len(dict["criteria"])):
#     for i in range(len(dict["subsystem"])):
#         subsystem = dict["subsystem"][i]
#         criteria = dict["criteria"][j]
#         #value_list = []
#         for k in range(len(weights)):
#             value = rating[i][j] * weights[k]
#             value_list[j][i][k] = value
#
#         #print(criteria, subsystem, value_list)
#         #print(subsystem, criteria, value)
#     print(value_list)
final_list = []
for j in range(len(dict["criteria"])):
    #final_list.append(dict["criteria"][j])

    for i in range(len(dict["subsystem"])):
        subsystem = dict["subsystem"][i]
        criteria = dict["criteria"][j]
        value_list = []
        for k in range(len(weights)):

            value = rating[i][j] * weights[k]
            value_list.append(value)

        final_list.append([dict["subsystem"][i], value_list])
        #print(criteria, subsystem, value_list)
            #print(subsystem, criteria, value)

print(final_list)

new_array = [final_list[0][1][8], final_list[1][1][8],final_list[2][1][8],final_list[3][1][8]]
# value = rating[0][0] * weights[6]
# print(value)
print(new_array)
original_value_list[0] = new_array
print(original_value_list)
new_sum = np.sum(original_value_list,axis=0)
#ranking_original = np.sort(original_value_list[0])
ranking_original = np.flip(np.argsort(original_sum))
ranking_new = np.flip(np.argsort(new_sum))

# print(ranking_original)
# print(ranking_original)
# print(original_value_list[0])

ranking_subsystem = []
for p in range(len(ranking_original)):
    choice = dict["subsystem"][ranking_original[p]]
    ranking_subsystem.append(choice)
print(ranking_subsystem)
print(new_sum)
print(original_sum)
new_ranking_subsystem = []
for p in range(len(ranking_new)):
    choice = dict["subsystem"][ranking_new[p]]
    new_ranking_subsystem.append(choice)
print(new_ranking_subsystem)


#print(dict["subsystem"][ranking_original[i]])

