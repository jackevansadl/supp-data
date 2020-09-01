for i in *; do sed -i '2s/.*/NumberOfCycles                1000/' ${i}/simulation_movies.input; done

for i in *; do sed -i '3s/.*/NumberOfInitializationCycles  1000/' ${i}/simulation_movies.input; done

for i in *; do sed -i '5s/.*/RestartFile                   yes/' ${i}/simulation_movies.input; done

for i in *; do sed -i '6s/.*/Movies yes/' ${i}/simulation_movies.input; done

for i in *; do sed -i '7s/.*/WriteMoviesEvery 10/' ${i}/simulation_movies.input; done

