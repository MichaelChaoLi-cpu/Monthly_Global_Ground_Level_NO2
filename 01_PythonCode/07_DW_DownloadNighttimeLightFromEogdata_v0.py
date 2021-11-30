# -*- coding: utf-8 -*-
"""
Spider Download Data from https://eogdata.mines.edu/nighttime_light/monthly/v10/2021/202106/vcmcfg/

Data: Monthly Nighttime Light 

Created on Tue Nov 30 11:07:40 2021

@author: M.L.
"""

from selenium.webdriver.support.ui import Select
import time
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
import pandas as pd

driver = webdriver.Chrome(ChromeDriverManager(version="95.0.4638.69").install())

year = 2021
month = ['01', '02', '03', '04', '05', '06',
         '07', '08', '09', '10', '11', '12']
i = 0
while i < 7:
    locationService = 'https://eogdata.mines.edu/nighttime_light/monthly/v10/' + str(year) + '/' + str(year) + month[i] + '/vcmcfg/'
    
    driver.get(locationService)
    driver.find_element_by_xpath(r'//*[@id="indexlist"]/tbody/tr[2]/td[2]/a').click()
    time.sleep(10)
    driver.get(locationService)
    driver.find_element_by_xpath(r'//*[@id="indexlist"]/tbody/tr[3]/td[2]/a').click()
    time.sleep(10)
    driver.get(locationService)
    driver.find_element_by_xpath(r'//*[@id="indexlist"]/tbody/tr[4]/td[2]/a').click()
    time.sleep(10)
    driver.get(locationService)
    driver.find_element_by_xpath(r'//*[@id="indexlist"]/tbody/tr[5]/td[2]/a').click()
    time.sleep(10)
    driver.get(locationService)
    driver.find_element_by_xpath(r'//*[@id="indexlist"]/tbody/tr[6]/td[2]/a').click()
    time.sleep(10)
    driver.get(locationService)
    driver.find_element_by_xpath(r'//*[@id="indexlist"]/tbody/tr[7]/td[2]/a').click()
    time.sleep(10)
    i = i + 1