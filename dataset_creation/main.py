import mechanism_generator
import statistic_module
import multiprocessing
import reaction_processor
import json


def mechanism_database_creation(num_parallel, start_step, end_step, radius, debug):
    # list_of_func.database_construct(num_parallel)
    # Your database should go here
    mechanism_generator.database_construct(r'', num_parallel, start_step, True)
    # IN CASE OF CREATING DATABASE


    with open(r'database.json') as db:
        mechanism_database = json.load(db)
        for num in range(num_parallel):
            with open(r'mech_db/mechanism_database_%d.json' % num, 'a') as write:
                json.dump(mechanism_database, write)

    if debug:
        reaction_processor.process_reaction(0, start_step, end_step, radius)
    else:
        # Create a list of arguments to pass to the worker processes
        args_list = [(num, start_step, end_step, radius)
                     for num in range(num_parallel)]

        # Create a pool of 32 worker processes
        with multiprocessing.Pool(processes=num_parallel) as pool:
            # Apply the worker function to each argument in parallel
            pool.starmap(reaction_processor.process_reaction, args_list)
            pool.close()
            pool.join()

    statistic_module.finalize_statistics(num_parallel)
    final_stat = statistic_module.calculate_total_statistics(num_parallel)
    statistic_module.write_finalized_statistics(final_stat, r'Data/Stat.txt')


if __name__ == '__main__':
    mechanism_database_creation(50, 1, 3, 3, False)



