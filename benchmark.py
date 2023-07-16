import pymolpro
import pandas as pd


def method_constructor(methods: list[str], counterpoise=False):
    return {
        method:
            [
                (openclosed + method.lower() if 'PAPT' not in method else openclosed + 'hf; papt; mp' + method[-1]) + (
                    '; counterpoise' if counterpoise else '')
                for openclosed in ['', 'U']]
        for method in methods
    }


def benchmark(project_name, database, methods, bases, max_electrons=None, open_shell=None, reference=None,
              backend='local',
              parallel=None, precision=2, units='kJ/mol', verbose=False, check=False, initial='symmetry,nosym'):
    answers = {}
    db = pymolpro.database.load(database).subset(max_electrons=max_electrons, open_shell=open_shell) if type(
        database) == str else database
    if verbose: print('size of database =', len(db))
    results = {}
    answers['bases'] = bases
    answers['units'] = units
    if units == 'kJ/mol': pd.options.display.float_format = '{:,.2f}'.format
    # print('answers[bases] set to',answers['bases'],type(answers['bases']))
    for method in methods:
        results[method] = {}
        for basis in bases:
            if verbose: print("run", method, basis)
            results[method][basis] = pymolpro.database.run(db, methods[method], basis, location=project_name,
                                                           initial=initial,
                                                           backend=backend, parallel=parallel, check=check)
            if check: print("completed run", method, basis)
    for method in methods:
        for k, v in results[method].items():
            for kp, vp in v.failed.items():
                print('failed:', method, k, kp, vp.filename())

    if check: print("after failure check")
    if len(results[list(methods.keys())[0]]) > 1:
        for method in methods:
            for resul in pymolpro.database.basis_extrapolate(results[method].values(), results['HF'].values()):
                results[method][resul.basis] = resul
            # print (results[method])
            for basis in results[method]:
                # print('look for basis ',basis,' in ',answers['bases'])
                if basis not in answers['bases']: answers['bases'].append(basis)

    answers['results'] = results

    pd.set_option('display.precision', precision)
    # print('reference',reference)
    reference_db = db if reference is None else results[reference][answers['bases'][-1]]
    answers['reference'] = reference_db
    # print('reference_db',reference_db)
    # print('results',results)
    # print( pymolpro.database.analyse([results[method][answers['bases'][-1]] for method in methods], reference_db, units))
    anal = pymolpro.database.analyse([results[method][answers['bases'][-1]] for method in methods], reference_db, units)
    for k, v in anal.items():
        # print(k,'energies' not in k,reference)
        # v.columns = pd.MultiIndex.from_arrays([methods.keys(), [c[1] for c in v.columns]])
        # if reference and 'energies' not in k: v.drop(columns=(reference,v.columns[-1][1]),inplace=True)
        if type(v) == pd.DataFrame:
            v.columns = methods.keys()
            if reference and 'energies' not in k: v.drop(columns=reference, inplace=True)
            answers[k] = v
    return answers


def all_project_directories(benchmarks):
    if type(benchmarks) != list: return all_project_directories([benchmarks])
    import os
    return [n.project_directory for bb in benchmarks for m in list(bb['results'].values()) for n in
            m.values() if n.project_directory is not None]


def unused_project_directories(benchmarks):
    import os
    dir = os.path.dirname(all_project_directories(benchmarks)[0])
    return [dir + '/' + x for x in os.listdir(dir) if
            x not in [os.path.basename(d) for d in all_project_directories(benchmarks)]]


def latex_table(df: pd.DataFrame, caption: str, reference_method=None, precision: int = 2, omitted_methods=[]) -> str:
    methods_pruned = [method for method in df.columns if method not in omitted_methods and method != reference_method]
    return f'{{\\ifx\\tablecaption\\undefined\\def\\tablecaption{{{caption}}}\\fi\n\\ifx\\toprule\\undefined\\def\\toprule{{\\hline\\hline}}\n\\def\\midrule{{\\hline}}\n\\def\\bottomrule{{\\hline\\hline}}\\fi % or \\usepackage{{booktabs}}\n' + \
        '\\newcolumntype{r}{>{\\raggedleft\\arraybackslash}p{3.0em}}' + \
        df[methods_pruned].style.format(
            precision=precision).to_latex(hrules=True, multicol_align='c', caption='\\tablecaption') + '}'


def violin_plot(bm, omitted_methods=['HF'], reference_method=None, zero_hline=True):
    ref_meth = list(bm['results'].keys())[-1] if reference_method is None else reference_method
    import matplotlib.pyplot as plt
    methods_pruned = [method for method in bm['results'] if method not in omitted_methods and method != ref_meth]
    fig, pane = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(6, 6))
    data = bm['reaction energy deviations'][methods_pruned].to_numpy()
    pane.violinplot(data, showmeans=True, showextrema=True, vert=True, bw_method='silverman')
    pane.set_xticks(range(1, len(methods_pruned) + 1), labels=methods_pruned, rotation=-90)
    pane.set_title(bm['bases'][-1])
    pane.set_ylabel('Error relative to ' + ref_meth + ' / ' + bm['units'])
    if zero_hline: plt.axhline(color='gray', linestyle=':', linewidth=0.5)
    plt.close()
    return fig


if __name__ == "__main__":
    bm = benchmark(
        'S66',
        'GMTKN55_S66',
        method_constructor(['HF', 'MP2', 'MP3', 'MP4', 'CCSD', 'PAPT2', 'PAPT3', 'PAPT4', 'CCSD(T)'],
                           counterpoise=True),
        # ['aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ'],
        ['aug-cc-pVDZ', 'aug-cc-pVTZ'],
        max_electrons=24,
        verbose=True,
        reference='CCSD(T)',
    )

    print(bm['reference'])

    print(bm['reaction energy deviations'])

    violin_plot(bm).show()
