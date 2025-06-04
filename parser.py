import re
from collections import defaultdict





def ParseStringsToDictsX(string_list):
    result = defaultdict(dict)
    for item in string_list:
        match = re.match(r"[A-Za-z]\[(\d+),(\d+),(\d+)\]", item)
        k_id, k_orig, k_dest = map(int, match.groups())
        if(not result[k_id]):
            result[k_id] = defaultdict(dict)
            max_a_id = -1
        else:
            max_a_id = max(result[k_id].keys())
        result[k_id][max_a_id+1] = (k_orig, k_dest)
    return list(result.values())


def ParseStringsToDictsY(string_list):
    result = defaultdict(dict)
    for item in string_list:
        match = re.match(r"[A-Za-z]\[(\d+),(\d+),(\d+),(\d+)\]", item)
        k_id, r_id, r_orig, r_dest = map(int, match.groups())
        if(not result[r_id]):
            result[r_id] = defaultdict(dict)
            max_a_id = -1
        else:
            max_a_id = max(result[r_id].keys())
        result[r_id][max_a_id+1] = (r_orig, r_dest)
    return list(result.values())

