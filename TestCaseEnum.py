from enum import Enum

class TEST(Enum):
    """
    Test case identification.
    """
    CASE_1 = 1 # paper case 1
    CASE_2 = 2 # paper case 2
    CASE_3 = 3  # jonas cross, order changes trains
    CASE_4 = 4 # jonas 3x3 k=1, r=2
    CASE_5 = 5 # jonas 3x3 k=2, r=2
    CASE_6 = 6 # jonas 3x3 k=2, r=2, diagonal arcs
    CASE_7 = 7 # jonas latvia simulation, no intermediate stops, r=3, k=3, same orig&dest for k and r
    CASE_8 = 8 # jonas latvia simulation, no intermediate stops, r=1, k=1, duplicate node on train route - ERROR
    CASE_9 = 9 # jonas latvia simulation, no intermediate stops, r=1, k=1, no duplicate node on train route
    CASE_10 = 10 # same as case_8, but train origin nodes doubled: n+1000=origin
    CASE_11 = 11 # same as 3, but 3 orders instead of 1
    CASE_12 = 12 # same as 3, but extended branches
    CASE_13 = 13 # latvia full
    CASE_14 = 14 # same as 8 but more nodes, k1 r4
    CASE_15 = 15 # same as 14 but more orders and trains
    CASE_16 = 16 # same as 8, k4, r4
    CASE_17 = 17 # same as 8, k1, r1, insufficient capacity
    CASE_18 = 18 # same as 8, k2, r2, specific transshipment terminals

class ROUTEPLOT(Enum):
    """
    Enum for plotting routes.
    """
    VEHICLE = 1
    ORDER = 2