def remove_percent_symbol(value):
    """
    Remove the '%' symbol from a percentage value.
    Args:
    - value (str): The percentage value with the '%' symbol.
    Returns:
    - value_without_percent (str): The percentage value without the '%' symbol.
    """
    return value.replace('%', '')

def test_answer() :
    yield1 = remove_percent_symbol("76%")
    assert yield1 == "76"

    yield2 = remove_percent_symbol("74.3%")
    assert yield2 == "74.3"

    yield3 = remove_percent_symbol(" ")
    assert yield3 == " "