# Retrieve hierarchy from wormmine
# This version retrieves just a single child
# Requires intermine installed: $ easy_install intermine
# cf. http://intermine.wormbase.org/tools/wormmine/query.do for query construction
#-------------------------------------------------------------------------------
# USAGE: python RetrieveHierarchy.py > hierarchy.csv
#-------------------------------------------------------------------------------

# Get intermine service
from intermine.webservice import Service
service = Service("http://intermine.wormbase.org/tools/wormmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("AnatomyTerm")

# Specify the output columns:
query.add_view("name","synonym","primaryIdentifier","children.name", \
    "children.primaryIdentifier","children.synonym")

# Specify a custom sort order?:
# query.add_sort_order("AnatomyTerm.name", "ASC")

#-------------------------------------------------------------------------------
# Just print names and IDs:
for row in query.rows():
    print '{0}|{1}|{2}|{3}'.format(row["name"],row["primaryIdentifier"], \
                        row["children.name"],row["children.primaryIdentifier"])

# Names, synonyms, and IDs:
# for row in query.rows():
#     print row["name"],"(",row["synonym"],"): ",row["primaryIdentifier"],",",row["children.name"], \
#         "(",row["children.synonym"],") ",row["children.primaryIdentifier"]
