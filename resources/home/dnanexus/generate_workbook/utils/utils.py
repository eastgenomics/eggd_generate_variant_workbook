from csv import Sniffer


def is_numeric(value:str) -> bool:
    """
    Returns true if given value is in some form numeric

    Parameters
    ----------
    value : str
        string to check if can be cast to numeric type

    Returns
    -------
    bool
        True if value can be numeric
    """
    try:
        if str(value)[0] == '0' and not str(value)[1] == '.':
            # catch things such as 01 which aren't real numbers
            return False
    except Exception:
        # could raise lots of errors (i.e IndexErrors from slicing)
        # lazily catch all and continue
        pass

    try:
        float(value)
        return True
    except ValueError:
        return False


def determine_delimeter(data) -> None:
    """
    Attempt to determine delimeter in a given string using csv.Sniffer,
    will default to tabs if it can't be determined

    Parameters
    ----------
    data : str
        data to check for delimeter

    Returns
    -------
    delimeter : str
        delimeter inferred from given data
    """
    try:
        delimeter = Sniffer().sniff(''.join(data)).delimiter
    except Exception as error:
        print(
            "Error in determing delimeter from given data. Will default "
            f"to using tabs.\n\nError: {error}\n\n"
        )
        delimeter = '\t'

    return delimeter
