from pyclustering.cluster import cluster_visualizer
from pyclustering.cluster.xmeans import xmeans
from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
from pyclustering.utils import read_sample
from pyclustering.samples.definitions import SIMPLE_SAMPLES
from enum import IntEnum

# Read sample 'simple3' from file.
#X = read_sample(SIMPLE_SAMPLES.SAMPLE_SIMPLE3)
#X = read_sample("16class.txt")
#X = read_sample("11class.txt")
X = read_sample("1class.txt")
print(len(X))
# Prepare initial centers - amount of initial centers defines amount of clusters from which X-Means will
# start analysis.
amount_initial_centers = 2
initial_centers = kmeans_plusplus_initializer(X, amount_initial_centers).initialize()
# Create instance of X-Means algorithm. The algorithm will start analysis from 2 clusters, the maximum
# number of clusters that can be allocated is 20.
xmeans_instance = xmeans(X, initial_centers, 20, False)

xmeans_instance.process()

# Extract clustering results: clusters and their centers
clusters = xmeans_instance.get_clusters()
centers = xmeans_instance.get_centers()
#print(clusters)
#print(centers)
#print(xmeans_instance.__bayesian_information_criterion(xmeans_instance.__clusters, xmeans_instance.__centers))
print(len(clusters), len(centers))
xmeans_instance.show_bic(clusters, centers);
xmeans_instance.show_mndl(clusters, centers);

# export dataset
f = open("sample_data.csv", "w")
for datapoint in X:
    element = 0
    for xy in datapoint:
        f.write(str(xy))
        if element == 0:
            f.write(",")
        element = element + 1
    f.write("\n")
f.close()

# export cluster index
f = open("cluster_indx.csv","w")
for cluster in clusters:
    element = 0
    maxElement = len(cluster)
    for indx in cluster:
        f.write(str(indx))
        if element < maxElement-1:
            f.write(",")
        element = element + 1
    f.write("\n")
f.close()

# export centroid coordinates
f = open("centroid_info.csv","w")
for center in centers:
    element = 0
    for xy in center:
        f.write(str(xy))
        if element == 0:
            f.write(",")
        element = element + 1
    f.write("\n")
f.close()

print("k = "+ str(len(centers)))
# Visualize clustering results
visualizer = cluster_visualizer()
visualizer.append_clusters(clusters, X)
visualizer.append_cluster(centers, None, marker='*', markersize=10)
visualizer.show()
