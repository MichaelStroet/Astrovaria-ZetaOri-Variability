
import os

def getPaths():
    """
    Finds the path toward several important directories.
    Returns a dictionary.
    """
    dirs = {}

    code = os.path.dirname(os.path.realpath(__file__))
    dirs["code"] = code

    main = os.path.dirname(code)
    dirs["main"] = main

    fastwind = os.path.join(main, "fastwind")
    dirs["fw"] = fastwind

    results = os.path.join(main, "results")
    dirs["results"] = results

    equivalent_width = os.path.join(results, "equivalent_widths")
    dirs["EW"] = equivalent_width

    zetaOri_obs = os.path.join(main, "data", "asci")
    dirs["zetaOri_obs"] = zetaOri_obs

    zetaOri_rv = os.path.join(main, "data", "asci-rv")
    dirs["zetaOri_rv"] = zetaOri_rv

    zetaOri_Annelotte = os.path.join(main, "data", "Annelotte")
    dirs["zetaOri_Annelotte"] = zetaOri_Annelotte

    return dirs

if __name__ == "__main__":
    dirs = getPaths()

    for dir in dirs:
        print(f"{dir}:\n{dirs[dir]}\n")
