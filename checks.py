"""functions for error prevention, e.g. checking data inputs, including:
    - check_H_dates
    - check_LC """


def check_H_dates(le_start, le_end):
    """check to make sure that end date is later than start date"""
    if le_start <= le_end:
        return
    else:
        print("The start date is after the end date! Check your dates in the Hydro workbook.")
        exit()


def check_LC(le_min, le_max):
    """checks that the minimum area available for land conservation is less than the max area available"""
    for i in range(len(le_min)):
        if le_min[i] > le_max[i]:
            print("Error: cannot have minimum area > maximum area for land conservation.")
            exit()

