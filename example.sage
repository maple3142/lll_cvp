from lll_cvp import *
from functools import partial, wraps

is_in_example = False


def example(fn):
    @wraps(fn)
    def wrapped(*args, **kwargs):
        global is_in_example
        if not is_in_example:
            print("=" * 40)
            print(f"Running {fn.__name__!r}")
            print(fn.__doc__)
            print("=" * 40)
        is_in_example = True
        fn(*args, **kwargs)
        is_in_example = False
        if not is_in_example:
            print()

    return wrapped


@example
def example1():
    """
    HITCON CTF 2019 Quals not so hard RSA
    """
    # copied from https://github.com/rkm0959/Inequality_Solving_with_CVP/blob/main/Example%20Challenge%204%20-%20HITCON%20CTF%202019%20Quals%20-%20not%20so%20hard%20RSA/solve_challenge_4.sage
    ## Example 4 : HITCON CTF 2019 Quals not so hard RSA

    ## d is 465 bits

    data = [
        (
            61608417975397048843788515638593839325111098880518441270527767841153782846066445099077365303960932518098100778959123136871039627996767023258612684873083420234538156646585282154245553305607644427220207313162116929585370583379703086997585339296145409828300576290109728682441066135201997295424597733433471586151,
            60032368056605168202792776655067640210910930719068898740685488293392455428589220656480049668823171895161714617099267690524276842795335016835073541061601545456195765907623303970146386500563913899580929779870429659650425339185233299118860275385880359287380867251468679962048998842668813298548390941601249105855,
            0x1E4433543AD3EAB1D5A5490E33EE98C34785945C7B69DD0FD0A371C28E5FF45F6627AD0559D9837FD6439367543FF5670F4DF4FD36CBEE75950DB62E51811F98E3F34DB66B07196A5DFBD9867952D8E6D67C43BECF086087181E5F78582E98945E5C8C08D754B998EF01E836729F9620CDCD2CC8AAE9CB4BF3D8E4BEEC3CA8FD,
        ),
        (
            52084595054768217522676979342755393306305099169414947960508049057119329537162079071100773540172780699842974838973453517862170568297372572652378982794735495836797191260523386985555538123219596180065060457770169661533302156160201909218943628670173938802399175928847179180982290051008255643478089280887936398779,
            15937444970326662770243305998311639639802081677519521519037892156989918335411343952028567001158781869582357356721336548356177945486289103616332672857562241745733313956089488246637981047367689313876169403573606972534488436195086998111247112950511965259280728899941655052549135516189815880662766290132607874671,
            0x15DBECED76AB710B84982E67A839846CFE38CFCCCE4DACC585DF0E38D695E1C84EAC7281D8C83B3C6EEBFAE7C27D91496B19DE120374D08CBCCC251A464C2F7BB10FB9D1C1F13B78F1FBFE0CE37D01350978DCB192C92DB43560A9CD81481AA2E2D41A8C5B3C67D3E3AE4BE50DDF37A5A0193DA6F4BEFE71D5348BD7820F1B0A,
        ),
        (
            121675354261226402523384817752122501670754379029920951567545832861118347020889141951034949457467042630795199347469665294677977729959566390358589470405088739543312116123816666440234661644847943610733959265939955886506334607329368744886743083353218982597246897307234102346047327550059519916639261619453107169923,
            1958149621008109700386193021020256359555810444022322757777842537494487439895477039590251808463946583928481366538370998263830913239177109491659797444270949612728937519312683847609162235917948298810271512469482956117395111528004250206351805086758526770795159863075378020975086922955262230095667701740508170151,
            0x5B445089B4578B115B0293CC1922F5FDB784701FD533EEC7EC9BD7FDAD995BAEFB051B9793FF3DADC24FF8D5B52C89F9565F65409C58506C7E79DD787F8E388F019497461FA3DB0EAC5284D398F9E6A1B81C59BA74677CC01C38DDD6461DF029E1179F3EE63CCCE20F3090835B7BD7B4DE25C38CCDDF2EC622CD41FAFB49D68F,
        ),
        (
            70820434096887624688036248070441718528107270792728727307197623356369639890276828895683705436125317302252834211775910545775001651922925694194466935425170969939539346172633115232030274859268960310846965871663926227518137713901150141826470260092611613639370170779217268337906342955347156686598981918675591938393,
            6264211827826908864953215815104077572906230669747226774603059704547415338199956437981145121700015020477626176620200121575740530313039227764537541481272011717745292622796843631960318813842515397395192806936521180850750841260559451018598913158548097304919304754876270597286947644956943158352022152801127867879,
            0x4F56516D8E197A8FD6ED76433A836635DD5AD1247BE9CCBB5A88EA940F5746221B4BB5B60EF925019A11D36FBE8F1D948E9AFB4305937D7E017148B8BA324682D60ED6FB7F3DE80031432023BBEF81A96D0C1BD26C3A728BC6FDEBBD48BF86B93325584C1A1386E26374C9747754B858A01B73996CEA8FE4EDD2A8130504E63E,
        ),
        (
            82396665631738285668995082087133930047188890932442133336046256534530991902128696466596278088102928917626464109384091264437455382690583003229076491363732222934817133536320704449557698925700841043105303694696203472531025908797520768019280257689521302427077849286757365739653728248123114986191011560571298134089,
            77541341309162852568860774868359137598587882474829780783306540898574011860779494674673722444562346968736216806881154962981695851598965376893425170303668270087989094501452417979429678869102952038600630142304872749698942475289176861137563927160985001506014089253274664444041660282426708388957971134224538326179,
            0x6B35DFFADD327C0999EFC909A3C1D1482C6A286808801095EEC5EA88224467881B4081C1AEF02C273CC5DC4D3505DC50FCF4CE60052B6A5A9DC005FAAD4709FBCCA254C6FC1C552D51C8E15FBD8FC404B0136758E1BA57F6C04B1049E303A43AE60C1EB0B671289F6689D1CD104C407549C1BCFC081A28A3F3324A1611E81555,
        ),
        (
            64885222129962661919689742957615146346855998523264787345799887210230045803423356080820154546050540728264000944692170729880161760807788699983366314269171423680987673788784023826885067996007133127166699735414920427997399513922777341273346196540175237091469555617281522541592407225697228219483666768990986901207,
            47625666889348051674457306461777570860477708163951567694477110559722255361844497439523737482859577232793119885759760321686098657292699849459852343686370029321619072315457562131610624443702007471950020751311592439673932509558946203816066868221203483167079506399876187240483994566182748627389022243608934595071,
            0x0F79419A86361D668CED50FBCFA521C5117F7B2F72C3CC248DBA4B2B1C9E7DC3A32FA500F5167A0360AEEFCB8DAA973BC67B0537D641617C2D96D98EE27179DE1644B480C120B2D14054DB032B42AF33B23F8182F72CC8E752D6F89D556800F637BC492BD8EB2FB294EEE61BC3C55677066B6962C4D2A1F1896703D446D903F3,
        ),
        (
            130307595686638523389042871138355252466828093625912424932382409748014813151912412889866698513547588105179569464359905871041351306187955332098310075195714887141180801001191398677311128761953690562878908997464309620384635922453765468745206262477334243233099872249671186622159533482223573972729351847574480258349,
            15517881429792452627724339825487038872077384939873021432131971858090290889497162855962736652640866521963165312583164013911162455865793567436833813352448640307756786561236246631146444590597902995591016928432988304077340670079503467394770901559567924809126517090243024890436892801283564005789256977822007581731,
            0x69A8E0BB49427F3465C9F3A41A22E047BD348C886CDB07264A321F8B890BB48CD7A878E43C1EB4B2F496DAFB50677B3EEA032C8F7F2EE59695398C56CC3B183BB6E2F1AB2D5D633461A6592CA0C98C0F2DAE4100D6AAAF0D1166BDD46BEB68B07D9C9CF6C5A92EF7DB019DC065B3A8300CA25B50CBA51A3294D954C178BF3770,
        ),
        (
            65386448419573832864151040666988491321490532636782967162806795276671017479566436411594017000387653366263378617189553041765586083217125727980734225153709370447270596363822524978426540069479632490507360491377612488876553715272119578498102854909764339676652952568021321118127618560269626549257296451954041925283,
            62417639911670877600600895776941278119707067464948560335387709799367122896567813566638216057171867961437006927011778910808189676129022776796728833211435603542364704330783289218396426543693678813151580911213239650850603907199690825764030369692366216968117035006819175280930722880390482701178418013436415892711,
            0x2231E71788FE0FEC78ADA4F24A58BBAA5B1A0FB48DC94B18BF7409637DBCFD8E2755D80B5E4B309E20A69474A6E2EB245E40C1B3F81C6151DBF7E6871E90A9616A32F4A8443D30BFAAFF5D711BE39BC8F38C13234AF9FD867F508B8A6B097A49FB07A1392EC08C38F0939620CF644EC6631F7566B6BC7A1F1D01A9C736ECEDA4,
        ),
        (
            131618812326147470215836136712095159211663684841466199972335887950664101678748935713203393931245541362772506324940879455504378545672259298894346211770983161798960234145317156224604949365110464634198678301644727985645228605590810190481288852412139519629077007437022950591320600450976104529132138396249413028099,
            14318453997383587408499979509282945054981178882287936277795320451723285671391694914929520724692059461777785480850364952738183846512279857119933447174333240691522012746191524101263192457454032837995702952840654708026147079321112397181981380295129696518148352311994008219805602114421675331461920408832158576071,
            0x71A370974F020D338071E2B0480F38F2D5B488252F5EB206636D6DCD3C93EA586A507C29D2E611A9A8D5D0F849B913116B37E69345D6EAE71CF87CDB6B74E75BCCDB374C372562777A727F07E4272A5FF17B451D582074B565879453D028C5E9B91CB67AE923F0F09492E508422E65C1DAADA0564D91A9FEEC094AF77AB190FA,
        ),
        (
            99041226655569332839951771030354960649881033648844611842019981984694036670092050113949713459076512100815903654289240773362239468428010238444408944822516190193861719380849620966232652503444079295871086550223067647300451354087943235559953597249075277760677450645744584993036921709483731953822197438848939991653,
            45917328654051873455626365238806652384056265616606200025887647202756379007355041189714541493038574645714670682101129613041399879280497291675656671383417575172090283625761438606371062142083212903854592596935358334738253893062615061401134834123551055465985828525720791523020773994634134135284036947641526947783,
            0x336E810E66EEAC255F1714B290187773DE92EC044B654A4CC0EEF2070E9D8D86FA37D424B2E4B71FD135E634034E06DBC31E77BD13D470AA5D3E5DCC5E689B565178A4A445D4560F5569967315A0BBE163EF83063429EC8D1FBBA1A8DBD9E30A2A67E17793F8B19E7B184C3B989DADC85B7F060E3DAEFB66FF803C98BDA8862F,
        ),
    ]

    ## each data has n, e for fixed d
    ## ed = k(n-p-q+1) + 1 -> ed + kn == k(-p-q+1) + 1
    ## construct a bound on k(-p-q+1) + 1
    ## 2 sqrt(n) <= p + q <= 3 sqrt(n/2)
    ## (e * 2^464 - 1) / (n - 2sqrt(n) + 1) <= (e * d - 1) / (n - 2sqrt(n) + 1) <= k
    ## k <= (e * d - 1) / (n - 3sqrt(n/2) + 1) <= (e * 2^465 - 1) / (n - 3sqrt(n/2) + 1)
    ## combine these to get a decent bound for k(-p-q+1) + 1

    ## 11 variables, d, and k for each 10 equations
    ## 11 equations, bound on d and each bound on ed + kn

    # build matrix
    M = matrix(ZZ, 11, 11)
    lb = [0] * 11
    ub = [0] * 11

    # encode d
    M[10, 10] = 1
    lb[10] = 2**464
    ub[10] = 2**465

    # encode ed + kn
    for i in range(0, 10):
        M[10, i] = data[i][1]  # e * d
        M[i, i] = data[i][0]  # k * n
        low_sum = int(2 * (data[i][0] ** 0.5))
        high_sum = int(3 * ((data[i][0] // 2) ** 0.5))
        low_k = (data[i][1] * (2**464) - 1) // (data[i][0] - low_sum + 1)
        high_k = (data[i][1] * (2**465) - 1) // (data[i][0] - high_sum + 1)
        lb[i] = high_k * (-high_sum + 1) + 1
        ub[i] = low_k * (-low_sum + 1) + 1

    res = solve_inequality(M, lb, ub)
    recovered_d = res[10]

    n = data[i][0]
    enc = data[i][2]

    ptxt = pow(enc, recovered_d, n)
    print((int)(ptxt).to_bytes(128, byteorder="big"))


@example
def example2():
    """
    ACSC Share The Flag
    """
    # modified from https://github.com/rkm0959/Inequality_Solving_with_CVP/blob/main/Example%20Challenge%207%20-%20ACSC%20Share%20The%20Flag/solve_challenge_7.py
    p = 251
    X = bytes.fromhex("02d4623be12c8f01cb2ebe5f837c1d")
    Y = bytes.fromhex("bbdc06ceb34da7b16336b007dc5492")
    X2 = bytes.fromhex("2fb9e753b237e68d35e266b0f01c9e")
    Y2 = bytes.fromhex("20c0be9140f5a33d71b9e82f8f9409")
    X3 = bytes.fromhex("f42e3ee10edeade0a3804a22e86a63")
    Y3 = bytes.fromhex("c7224da73d9d96254f94136d9a65f1")
    X4 = bytes.fromhex("37c9b07870283dd3f6198c46f027dd")
    Y4 = bytes.fromhex("8101a88a365526e8faf417b79599a0")
    X5 = bytes.fromhex("b0342cb7b3f5a022d927f9019a1bf3")
    Y5 = bytes.fromhex("e2666d892955494775aa3c96c441f5")
    X6 = bytes.fromhex("e56bf4f9e746252dbacb93a0a95087")
    Y6 = bytes.fromhex("cbb43831857333b2c4663ba2c9189a")
    X7 = bytes.fromhex("99ca36b1633cf3d903d8e6291f1bdc")
    Y7 = bytes.fromhex("25180068651818171d10422dbdb395")

    M = Matrix(GF(p), 105, 128)
    vec = []
    for i in range(105):
        x, y = 0, 0
        if i < 15:
            x = int(X[i])
            y = int(Y[i])
        elif i < 30:
            x = int(X2[i - 15])
            y = int(Y2[i - 15])
        elif i < 45:
            x = int(X3[i - 30])
            y = int(Y3[i - 30])
        elif i < 60:
            x = int(X4[i - 45])
            y = int(Y4[i - 45])
        elif i < 75:
            x = int(X5[i - 60])
            y = int(Y5[i - 60])
        elif i < 90:
            x = int(X6[i - 75])
            y = int(Y6[i - 75])
        elif i < 105:
            x = int(X6[i - 90])
            y = int(Y6[i - 90])

        vec.append(y)
        for j in range(16):
            M[i, j] = (x**j) % p
        if i < 15:
            for j in range(16):
                M[i, j + 16] = (x ** (j + 16)) % p
        elif i < 30:
            for j in range(16):
                M[i, j + 32] = (x ** (j + 16)) % p
        elif i < 45:
            for j in range(16):
                M[i, j + 48] = (x ** (j + 16)) % p
        elif i < 60:
            for j in range(16):
                M[i, j + 64] = (x ** (j + 16)) % p
        elif i < 75:
            for j in range(16):
                M[i, j + 80] = (x ** (j + 16)) % p
        elif i < 90:
            for j in range(16):
                M[i, j + 96] = (x ** (j + 16)) % p
        elif i < 105:
            for j in range(16):
                M[i, j + 112] = (x ** (j + 16)) % p

    vec = vector(GF(p), vec)

    bas = M.right_kernel().basis()
    print(len(bas))
    v = M.solve_right(vec)

    # v + bas -> all in 97 ~ 122

    M = Matrix(ZZ, 151, 151)
    lb = [0] * 151
    ub = [0] * 151

    for i in range(23):
        for j in range(128):
            M[i, j] = int(bas[i][j])
        M[i, 128 + i] = 1
    for i in range(128):
        M[23 + i, i] = p
    for i in range(128):
        if i >= 16:
            lb[i] = int(97 - int(v[i]))
            ub[i] = int(122 - int(v[i]))
        else:
            lb[i] = int(32 - int(v[i]))
            ub[i] = int(128 - int(v[i]))
    for i in range(23):
        lb[i + 128] = 0
        ub[i + 128] = p

    from functools import partial

    res = solve_inequality(
        M,
        lb,
        ub,
        cvp=partial(kannan_cvp, reduction=lambda M: M.BKZ(block_size=20), weight=251),
    )
    print(vector(lb))
    print(res)
    print(vector(ub))

    flag = ""

    for i in range(16):
        flag += chr((res[i] + int(v[i]) + 251 * 30) % 251)

    print("ACSC{" + flag + "}")
    # ACSC{wOAdvfst41xJzG6r}


@example
def example3():
    """
    PBCTF 2021 - seed me
    """
    # modified from https://blog.maple3142.net/2021/10/11/pbctf-2021-writeups/#seed-me
    from operator import xor

    class JavaRNG:
        # about the detail of java 17 rng
        # https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/util/Random.html
        def __init__(self, seed):
            self.seed = seed

        def next(self):
            self.seed = self.seed * 0x5DEECE66D + 0xB
            return self.seed

    Z = Zmod(2**48)
    P = PolynomialRing(Z, "s")
    s = P.gen()

    aa = []
    bb = []
    zz = []
    rng = JavaRNG(s)
    for _ in range(16):
        for _ in range(2047):
            rng.next()
        z = rng.next()
        # print(z)
        zz.append(z)
        b, a = z.change_ring(ZZ)
        aa.append(a)
        bb.append(b)
        # print(((ZZ(z(xs)) >> 24) / (1 << 24)).n())

    M = 2**48
    B = block_matrix(
        [[matrix([1]), matrix(aa)], [matrix(len(aa), 1), matrix.identity(len(aa)) * M]]
    )

    # manually changing parameters...
    vlb = [M - 3100000000000 for _ in range(len(bb))]
    vlb[0] = M - 2900000000000
    vub = [M - 2100000000000 for _ in range(len(bb))]
    vub[-1] = M - 2400000000000

    lb = [0] + [v - b for v, b in zip(vlb, bb)]
    ub = [2**48] + [v - b for v, b in zip(vub, bb)]
    res = solve_inequality(matrix(B), list(lb), list(ub))

    s = ZZ(res[0])
    for z in zz:
        r = ZZ(z(s))
        o = ((r >> 24) / (1 << 24)).n()
        print(o, o > 7.331 * 0.1337)

    print(
        xor(s, 0x5DEECE66D)
    )  # Java RNG will xor your seed with 0x5DEECE66D when setting seed
    # known good seed: 272404351039795


@example
def example4():
    """
    Simple example of underconstrained equations
    """
    n = 10
    pub1 = random_vector(ZZ, n, x=1, y=2**256)
    pub2 = random_vector(ZZ, n, x=1, y=2**256)
    secret = random_vector(ZZ, n, x=1, y=2**64)
    t1 = pub1 * secret
    t2 = pub2 * secret
    print(
        solve_underconstrained_equations(
            matrix([pub1, pub2]).T, vector([t1, t2]), [0] * n, [2**64] * n
        )
    )
    print(secret)


@example
def example5(cvp=kannan_cvp):
    """
    ImaginaryCTF Round 30 - Easy DSA: LCG
    """
    # https://github.com/maple3142/My-CTF-Challenges/tree/master/ImaginaryCTF/Round%2030/Easy%20DSA:%20LCG
    from fastecdsa.curve import secp256k1
    from hashlib import sha256
    from Crypto.Cipher import AES

    G = secp256k1.G
    q = secp256k1.q
    # fmt: off
    p = 9927040122486684509203958106419420141058188722199373989012953585197167125223276141324574147521754273735827724127795605194092299982453828901469369136978219
    sigs = [(98078224267884884220741740422077019843954009281647502734600509731511013529371, 54523865988310606978987830048871561792183822750263202533230451076893555969316), (104372973739209868434840748268723094332969140159620819033951611727659419363988, 39660851627725578124743718742328950528148285144862142963822549722002689280409), (103919709879086855178251181244489637133481828592253195107866903154222896468253, 35031204282583023574328215246485186335362731664384171126097342931654133207246), (63175283280752608661708773461972110889312169792285211062806717970617630555061, 34712080692439206749112321272818736084925608248138548106200594874651099131535)]
    ct = b'\xe6\x9c\xcaZ\x01\x90-\xa0\xbc8\xeb\xe4\xc6\xc7b\x16\xb9t++@\xc0\x0ce\t\x9e\xb5\x07p\xe49*\xb8\xce\xfe@\xea%\xc9\xd6\xefF\xf8\x7fQ\x9bg\xbd\x7f\xcf{h\\^\x11\xf9\xf5\xe8\x7f}\x94\xd3+\x06\x19.`\x84\x8d)\x1e\xdey\xe4 [\x9e'
    nonce = b'Z\x1c\xba\xbc\x95\\\xe1u'
    # fmt: on
    msgs = [
        b"https://www.youtube.com/watch?v=S8MJvhgjXBY",
        b"https://www.youtube.com/watch?v=wSTbdqo-j74",
        b"https://www.youtube.com/watch?v=dkYHgxfQZBA",
        b"https://www.youtube.com/watch?v=p8ET-m6y6VU",
    ]

    ss = []
    for m, (r, s) in zip(msgs, sigs):
        z = int.from_bytes(sha256(m).digest(), "big") % q
        ss.append((z, r, s))
    a, b = G.x, G.y
    syms = "d," + ",".join([f"k{i}" for i in range(len(ss))])
    R = ZZ[syms]
    d, *ks = R.gens()
    # collect equations
    eq_p = []
    eq_q = []
    for (z, r, s), k in zip(ss, ks):
        eq_q.append(s * k - z - r * d)
    eq_q = [f.resultant(g, d) for f, g in zip(eq_q, eq_q[1:])]
    for k, kk in zip(ks, ks[1:]):
        eq_p.append(a * k + b - kk)

    # solve!
    eqs = eq_p + eq_q
    mods = [p] * len(eq_p) + [q] * len(eq_q)
    lb = [0] * len(ks)
    ub = [2**512] * len(ks)
    ks = solve_multi_modulo_equations(eqs, mods, lb, ub, cvp=cvp)

    z, r, s = ss[0]
    k = ks[0]
    d = (s * k - z) * pow(r, -1, q) % q

    # get flag
    key = sha256(str(d).encode()).digest()[:16]
    cipher = AES.new(key, AES.MODE_CTR, nonce=nonce)
    print(cipher.decrypt(ct))


@example
def example6():
    """
    Simple example of underconstrained equations (non-linear equations)
    """
    n = ZZ(getrandbits(2048))
    roots = [ZZ(getrandbits(128)) for _ in range(3)]
    x, y, z = PolynomialRing(ZZ, ["x", "y", "z"]).gens()
    f = randrange(1, n) * x * y + randrange(1, n) * y * z + randrange(1, n) * z * x
    f -= f(*roots)
    f %= n
    g = randrange(1, n) * x**2 + randrange(1, n) * y**2 + randrange(1, n) * z**2
    g -= g(*roots)
    g %= n
    eqs = [f, g]
    bounds = {
        x: 2**128,
        y: 2**128,
        z: 2**128,
    }
    for solve in (
        solve_underconstrained_equations_general,
        solve_underconstrained_equations_general_v2,
    ):
        print(f"solving with {solve.__name__}")
        for monos, sol in solve(n, eqs, bounds):
            print(monos, sol)
            if sol[-1] < 0:
                sol = -sol
            if sol[-1] == 1:
                polys = [f.change_ring(QQ) for f in sol - monos if f]
                I = ideal(polys)
                print(I.variety())
                print(roots)


@example
def example7():
    """
    TSJ CTF 2022 - Signature
    """
    # https://github.com/maple3142/My-CTF-Challenges/blob/master/TSJ%20CTF%202022/Signature/README.md
    from fastecdsa.curve import secp256k1 as CURVE
    from Crypto.Cipher import AES
    from hashlib import sha256
    from operator import xor

    sigs = [
        (
            68628903551760154000300039815420920736833607628328728397865761689137575263983,
            15443603266495080692405049373908575802334630966756172054216234269624270469394,
            22981409647080659483954459639297160899607455746832364780169145673684744217331,
        ),
        (
            78718753887900409677937247756515988306958319672631190519685073607080826418294,
            54276644149626679124071775426283215897696024493497440097979013953212408253734,
            55041549251920287655053310452444903772819013365729781587201778245299568302064,
        ),
        (
            76095801129518609512809073275431964691465208338605553566026645469843430112342,
            13243268453692556572476611015151647125720859921369299900977957389460914346512,
            36340755739784353369950946724334315371788689739091478499529319569139128625422,
        ),
        (
            14730727965324371353827007106163274690512416783313793290093026161256262981745,
            61155679185854480315452772506628842842531730179984684102346013415628266033721,
            23105235892108839369406685343889200404631332972108556923344303689260664287078,
        ),
        (
            1873902112306026140517260545900147482389676444908353863350522632481552764819,
            7259111370509779279520618531388517380717524647694957730980169106713190467749,
            47964977855256430218836508671575141052659514288194710563910728384541064927004,
        ),
        (
            82140969634605912978034198684963146226227701184998942661897478141366806021323,
            83925169940944825117018170424384235795850122118608830570830266120233738605795,
            18779900738289102423955738918979596192495903397022523360371032408971788920887,
        ),
    ]
    msgct = b"\xd8n\x81\xcfrN=\xa8\xda\xa0pIR\xdb8\x98\xdb.s;\xecOG\xe2\x07\x99\x93\x1ah\x08G\xb2\x82\xe5\xef2\xe9\x92\x8b\xfaA\x8f\xa3\xa5\xcc8\x90\x95\xff?\xef\x1c\xfd\xffL(\xb46H\x0e/z/\xfe\x9eh[\xb3\xc6Y\x0e\x90\xe5\xdaU\xcc\x84\xed\x81\xcb\x1d<`]\xd2E\xd9\xa3\xa5u\xa3/b\xf3qz\x19\xfa\x8e$\xa3S\xac{,\xbcYIZ\xa9Z\xd6M-\x06Io>\xb2\xa35\x97=\x10]\xb36\xf6\xb2.\x8c\xd48\xe4"

    def recover_d(sigs):
        P = PolynomialRing(Zmod(CURVE.q), "d", 256)
        ds = P.gens()
        dd = sum([2**i * di for i, di in enumerate(ds)])
        polys = []
        for z, r, s in sigs:
            d_and_z = sum(
                [2**i * ((z & (1 << i)) >> i) * di for i, di in enumerate(ds)]
            )
            # fact: (a xor b) = a + b - 2 * (a and b)
            k = dd + z - 2 * d_and_z
            polys.append((s * k) - (z + r * dd))
        bounds = {d: 1 for d in ds}
        for monos, sol in solve_underconstrained_equations_general(
            CURVE.q, polys, bounds
        ):
            if abs(sol[-1]) == 1:
                print(monos)
                print(sol)
                dbits = (sol * sol[-1])[:-1]
                d = int("".join(map(str, dbits[::-1])), 2)
                return d

    d = recover_d(sigs)
    print(f"{d = }")
    key = sha256(str(d).encode()).digest()[:16]
    nonce = AES.new(key, AES.MODE_ECB).decrypt(
        bytes(map(xor, b"Congrats! This is your flag: "[:16], msgct[:16]))
    )[:8]
    cipher = AES.new(key, AES.MODE_CTR, nonce=nonce)
    print(cipher.decrypt(msgct))


@example
def example8():
    """
    SEETF 2023 - onelinecrypto
    assert __import__('re').fulmatch(r'SEE{\w{23}}',flag:=input()) and not int.from_bytes(flag.encode(),'big')%13**37
    """
    from Crypto.Util.number import bytes_to_long
    import re

    N = 23
    C = int.from_bytes(b"SEE{" + b"\x00" * N + b"}", "big")
    L = matrix.zero(N + 1, N + 1)
    L[0, 0] = 13**37
    for i in range(1, N + 1):
        L[i, 0] = 256**i
        L[i, i] = 1

    lb = vector([-C] + [ord("0")] * N)
    ub = vector([-C] + [ord("z")] * N)
    sols, basis = solve_inequality_ex(L, lb, ub)
    s0 = sols[0]
    not_flag = b"SEE{" + bytes(s0[1:][::-1]) + b"}"
    print(s0, not_flag)
    assert bytes_to_long(not_flag) % (13**37) == 0

    for sol in enum_ilp(sols[0], basis, lb, ub):
        val = bytes(sol[1:][::-1])
        flag = b"SEE{" + val + b"}"
        print(flag)
        assert bytes_to_long(flag) % (13**37) == 0
        try:
            if re.fullmatch(r"SEE{\w{23}}", flag.decode()):
                print("FOUND")
                break
        except UnicodeDecodeError:
            pass


@example
def example9():
    """
    Hidden Subset Sum Problem (https://eprint.iacr.org/2020/461.pdf)
    """

    p = 2**255 - 19
    n, m = 16, 64

    def gen_instance(p, n, m):
        alpha = random_vector(GF(p), n)
        X = random_matrix(GF(2), m, n).change_ring(ZZ).change_ring(GF(p))
        h = X * alpha
        return h, X, alpha

    h, X, alpha = gen_instance(p, n, m)

    # orthogonal lattice attack
    ortho = find_ortho(p, h)
    ortho2 = find_ortho(None, *ortho[: m - n])
    print(ortho2)
    assert ortho2.dimensions() == (n, m)
    # the entries of ortho2 are -1,0,1, but we want 0,1
    # so we can apply ilp based enumeration to solve it
    lb = vector([0] * m)
    ub = vector([1] * m)
    collect = []
    for sol in enum_ilp(None, ortho2, lb, ub):
        if sol == 0:
            continue
        collect.append(sol)
        if len(collect) >= ortho2.nrows():
            break
    X2 = matrix(collect).T  # X2.T and X.T are equivalent up to rows being permuted
    print(X2.T)
    assert X.T.is_permutation_of(X2.T)
    # so we can't guarantee the order the recovered alpha, but its sum is correct
    assert set(X2.solve_right(h)) == set(alpha)


@example
def example10():
    """
    LWE primal attack
    """
    q = 2**255 - 19
    n = 32
    m = 64
    e_bound = 10

    A = random_matrix(GF(q), m, n)
    s = random_vector(GF(q), n)
    e = vector(GF(q), [randint(-e_bound, e_bound) for _ in range(m)])
    b = (A * s + e).change_ring(ZZ)

    # LWE primal attack
    ATR = reduce_mod_p(A.T, q)
    s_rec = A.solve_right(kannan_cvp(ATR, b))
    print(s_rec)
    assert s_rec == s

    # Reduce it later also works
    ATR = reduce_mod_p(A.T, q, reduction=lambda x: x)
    s_rec = A.solve_right(kannan_cvp(ATR, b))
    print(s_rec)
    assert s_rec == s


@example
def example11():
    """
    ECDSA biased nonce
    """
    from fastecdsa.curve import secp256k1
    from fastecdsa.keys import gen_keypair
    from secrets import randbelow

    q = secp256k1.q
    G = secp256k1.G

    def sign(d, z, k):
        r = (k * G).x
        s = (z + r * d) * pow(k, -1, q) % q
        return r, s

    d, Y = gen_keypair(secp256k1)
    zs = [1, 2, 3]
    B = 0x13371337133713371337133713371337 << 128
    ks = [B + randbelow(1 << 128) for _ in zs]
    sigs = [sign(d, z, int(k)) for z, k in zip(zs, ks)]
    rs, ss = zip(*sigs)

    def ecdsa_eq(d, z, r, s, k):
        return s * k - z - r * d

    ds, k0, k1, k2 = polygens(Zmod(q), ("d", "k0", "k1", "k2"))
    eq0 = ecdsa_eq(ds, zs[0], rs[0], ss[0], k0)
    eq1 = ecdsa_eq(ds, zs[1], rs[1], ss[1], k1)
    eq2 = ecdsa_eq(ds, zs[2], rs[2], ss[2], k2)
    eqs = [eq0, eq1, eq2]
    sb = (B, B + 2**128)
    bounds = {ds: q, k0: sb, k1: sb, k2: sb}
    for monos, sol in solve_underconstrained_equations_general_v2(q, eqs, bounds):
        sol_dict = dict(zip(monos, sol))
        assert sol_dict[ds] % q == d
        print("found", sol_dict[ds] % q)


@example
def example_flatter():
    """
    Example 5 but explicitly using flatter reduction
    """
    example5(cvp=partial(kannan_cvp, reduction=flatter))


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.DEBUG)
    example1()
    example2()
    example3()
    example4()
    example5()
    example6()
    example7()
    example8()
    example9()
    example10()
    example11()
    example_flatter()
