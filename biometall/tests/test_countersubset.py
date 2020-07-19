import pytest
from biometall.modules import box

testdata = [
    ([], [], True),
    ([], ['HIS'], True),
    (['HIS'], [], False),
    (['HIS'], ['HIS'], True),
    (['HIS'], ['HIS', 'GLU'], True),
    (['HIS', 'GLU'], ['HIS'], False),
    (['HIS', 'GLU'], ['HIS', 'GLU', 'TYR'], True),
    (['HIS', 'GLU', 'TYR'], ['HIS', 'GLU'], False),
    (['HIS', 'GLU', 'TYR'], ['HIS', 'GLU', 'TYR'], True),
    (['HIS', 'GLU', 'GLU', 'GLU', 'TYR'], ['GLU', 'GLU', 'TYR', 'HIS', 'GLU'], True),
    (['HIS', 'GLU', 'GLU', 'GLU', 'TYR'], ['GLU', 'GLU', 'TYR', 'HIS', 'GLU', 'ASP'], True),
    (['HIS', 'GLU', 'GLU', 'ASP', 'GLU', 'TYR'], ['GLU', 'GLU', 'TYR', 'HIS', 'GLU'], False),
    (['HIS', 'GLU', 'GLU', 'GLU', 'TYR', 'HIS'], ['GLU', 'GLU', 'TYR', 'HIS', 'GLU', 'ASP', 'HIS'], True),
]
@pytest.mark.parametrize("list1,list2,result", testdata)

def test_countersubset(list1, list2, result):
    is_subset = box._counterSubset(list1, list2)
    assert is_subset == result
