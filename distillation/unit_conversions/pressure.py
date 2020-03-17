def psia_to_Paa(P):
    """

    :param P: pressure in psia
    :return: pressure in Pa (absolute)
    """
    return P / 14.5038 * 1e5


def Paa_to_psia(P):
    """

    :param P: pressure in Pa (absolute)
    :return: pressure in psia
    """
    return P / 1e5 * 14.5038