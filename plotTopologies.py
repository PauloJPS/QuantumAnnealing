### Importing Libraries 
import matplotlib.pyplot as plt
import dwave_networkx as dnx
import networkx as nx

### Creating topologies 
G_chimera = dnx.chimera_graph(3, 3, 4)
G_pegasus = dnx.pegasus_graph(3)
pos_chimera = dnx.chimera_layout(G_chimera)
pos_pegasus = dnx.pegasus_layout(G_pegasus)


### Plotting 
fig, ax = plt.subplots(figsize=(6,6))
nx.draw_networkx_nodes(G_chimera, pos_chimera, node_color='#377eb8', edgecolors='black', node_size=100, ax=ax)
nx.draw_networkx_edges(G_chimera, pos_chimera, ax=ax)
plt.axis('off')
ax.get_xaxis().set_visible(False) # this removes the ticks and numbers for x axis
ax.get_yaxis().set_visible(False) # this removes the ticks and numbers for y axis
plt.savefig('chimera.pdf', bbox_inches='tight',pad_inches = 0, dpi = 200)
plt.close()

fig, ax = plt.subplots(figsize=(6,6))
nx.draw_networkx_nodes(G_pegasus, pos_pegasus, node_color='#e41a1c', edgecolors='black', node_size=100, ax=ax)
nx.draw_networkx_edges(G_pegasus, pos_pegasus, ax=ax)
plt.axis('off')
ax.get_xaxis().set_visible(False) # this removes the ticks and numbers for x axis
ax.get_yaxis().set_visible(False) # this removes the ticks and numbers for y axis
plt.savefig('pegasus.pdf', bbox_inches='tight',pad_inches = 0, dpi = 200)



