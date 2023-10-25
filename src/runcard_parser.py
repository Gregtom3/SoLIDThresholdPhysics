# Program for parsing the user inputted runcard

def parse_runcard(args):
    runcard_path = args.runcard
    parameters = {}
    
    with open(runcard_path, 'r') as file:
        for line in file:
            # skip comments or empty lines
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            key, value = line.split(":")
            parameters[key.strip()] = value.strip()

    try:
        parameters["batch"] = args.batch
    except:
        parameters["batch"] = 0
        
    return parameters