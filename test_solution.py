from satellite import construct_satellite
from satellite import JD_START
from main import START_POS, START_VEL, TARGET_POINT, CALENDAR_TIME_START, convert_geo2eci 


def test_answer():
    """
    тест сравнивает координаты и скорости спутника, определенные алгоритмом  
    с заданными по условиям задания значениями в начальный момент времени
    """
    target_point = convert_geo2eci(TARGET_POINT, CALENDAR_TIME_START)
    
    params_time_begin = [START_POS['x'], START_POS['y'],
        START_POS['z'], START_VEL['vel_x'], START_VEL['vel_y'], START_VEL['vel_z']]
    time = 0  
    satellite = construct_satellite(START_POS, START_VEL, time, target_point)

    satellite.update_params(time)
    params = satellite.get_params()[1:7]
    print(params)
    print(params_time_begin)
    compare_flag = True
    for i in range(len(params_time_begin)):
        if not int(params_time_begin[i]) == int(params[i]):
            compare_flag = False

    assert compare_flag
