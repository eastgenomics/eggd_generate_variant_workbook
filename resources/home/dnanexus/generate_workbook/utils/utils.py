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
    return str(value).lstrip('-').replace('.', '').replace(
        'e-', '', 1).replace('e', '').isdigit()