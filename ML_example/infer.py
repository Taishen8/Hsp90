from torch.utils.data import TensorDataset, DataLoader
import math
import pickle 
import torch
import torch.nn as nn
from time import sleep
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu' )

class GlobalMaxPool1d(nn.Module):
    def __init__(self):
        super(GlobalMaxPool1d,self).__init__()
    def forward(self,x):
        output, _ = torch.max(x,1)
        return output

class ConvNN(nn.Module):
    def __init__(self,in_dim,c_dim,kernel_size):
        super(ConvNN,self).__init__()
        self.convs = nn.Sequential(
            nn.Conv1d(in_channels=in_dim, out_channels= c_dim, kernel_size=kernel_size,padding='same'),
            nn.ReLU(),
            nn.Conv1d(in_channels=c_dim, out_channels= c_dim*2, kernel_size=kernel_size,padding='same'),
            nn.ReLU(),
            nn.Conv1d(in_channels=c_dim*2, out_channels= c_dim*3, kernel_size=kernel_size,padding='same'),
            nn.ReLU(),
            #GlobalMaxPool1d() # 192
            )
    def forward(self,x):
        x = self.convs(x)
        return x

class Self_Attention(nn.Module):
    # input : batch_size * seq_len * input_dim
    # q : batch_size * input_dim * dim_k
    # k : batch_size * input_dim * dim_k
    # v : batch_size * input_dim * dim_v
    def __init__(self,input_dim,dim_k,dim_v):
        super(Self_Attention,self).__init__()
        self.q = nn.Linear(input_dim,dim_k)
        self.k = nn.Linear(input_dim,dim_k)
        self.v = nn.Linear(input_dim,dim_v)
        self._norm_fact = 1 / math.sqrt(dim_k)
        
    
    def forward(self,x):
        Q = self.q(x) # Q: batch_size * seq_len * dim_k
        K = self.k(x) # K: batch_size * seq_len * dim_k
        V = self.v(x) # V: batch_size * seq_len * dim_v
         
        atten = nn.Softmax(dim=-1)(torch.bmm(Q,K.permute(0,2,1))) * self._norm_fact # Q * K.T() # batch_size * seq_len * seq_len
        
        output = torch.bmm(atten,V) # Q * K.T() * V # batch_size * seq_len * dim_v
        
        return output

class CAMP(nn.Module):
    def __init__(self):
        super(CAMP,self).__init__()
        #self.config = config
        self.embed_seq = nn.Embedding(65+1, 128) # padding_idx=0, vocab_size = 65/25, embedding_size=128
        self.embed_ss = nn.Embedding(75+1,128)
        self.embed_two = nn.Embedding(7+1,128)
        self.pep_convs = ConvNN(512,64,7)
        self.prot_convs = ConvNN(512,64,8)
        self.pep_fc = nn.Linear(3,128)    
        self.prot_fc = nn.Linear(23,128)
        self.global_max_pooling = GlobalMaxPool1d()
        #self.dnns = DNN(config.in_dim,config.d_dim1,config.d_dim2,config.dropout)
        self.dnns = nn.Sequential(
            nn.Linear(640,1024),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(1024,1024),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(1024,512))
        
        self.att = Self_Attention(128,128,128)
        #c_dim
        self.output = nn.Linear(512,1)

    #@torchsnooper.snoop()
    def forward(self, x_pep,x_prot,x_pep_ss,x_prot_ss,x_pep_2,x_prot_2,x_pep_dense,x_prot_dense):

        pep_seq_emb = self.embed_seq(x_pep.long())#.type(torch.LongTensor))
        prot_seq_emb = self.embed_seq(x_prot.long())#.type(torch.LongTensor))
        pep_ss_emb = self.embed_ss(x_pep_ss.long())#type(torch.LongTensor))
        prot_ss_emb = self.embed_ss(x_prot_ss.long())
        pep_2_emb = self.embed_two(x_pep_2.long())
        prot_2_emb = self.embed_two(x_prot_2.long())
        pep_dense = self.pep_fc(x_pep_dense)
        prot_dense = self.prot_fc(x_prot_dense)
        

        encode_peptide = torch.cat([pep_seq_emb, pep_ss_emb, pep_2_emb, pep_dense],dim=-1)
        encode_protein = torch.cat([prot_seq_emb, prot_ss_emb, prot_2_emb, prot_dense],dim=-1)

        encode_peptide = encode_peptide.permute(0,2,1)
        encode_protein = encode_protein.permute(0,2,1)

        encode_peptide = self.pep_convs(encode_peptide)
        encode_peptide = encode_peptide.permute(0,2,1)
        encode_peptide_global = self.global_max_pooling(encode_peptide)

        encode_protein = self.prot_convs(encode_protein)
        encode_protein = encode_protein.permute(0,2,1)
        encode_protein_global = self.global_max_pooling(encode_protein)
        
        # self-attention
        pep_seq_att = self.embed_seq(x_pep.long())
        peptide_att = self.att(pep_seq_att)
        peptide_att = self.global_max_pooling(peptide_att)
        
        prot_seq_att = self.embed_seq(x_prot.long())
        protein_att = self.att(prot_seq_att)
        protein_att = self.global_max_pooling(protein_att)

        encode_interaction = torch.cat([encode_peptide_global,encode_protein_global,peptide_att,protein_att],axis=-1)
        encode_interaction = self.dnns(encode_interaction)
        predictions = torch.sigmoid(self.output(encode_interaction))

        return predictions.squeeze(dim=1)

ft_idx_lst = [50,639,689,1278,1328,1917,2067,15614]

amino_acid_dict = { "A": 1, "C": 2, "B": 3, "E": 4, "D": 5, "G": 6, "F": 7,
        "I": 8, "H": 9, "K": 10, "M": 11, "L": 12, "O": 13, "N": 14, "Q": 15,
        "P": 16, "S": 17, "R": 18, "U": 19, "T": 20, "W": 21, "V": 22, "Y": 23,
        "X": 24, "Z": 25 }
rev_map = {}
for k,v in amino_acid_dict.items():
    rev_map[v] = k

def infer(model,dataloader,device,info):
    model.eval()
    preds, seqs = [], []
    with torch.no_grad():
        avg_loss = 0
        for batch, X in enumerate(dataloader):
            X = X[0].to(device)
            X_pep = X[:,:ft_idx_lst[0]]
            X_prot = X[:,ft_idx_lst[0]:ft_idx_lst[1]]
            X_pep_ss = X[:,ft_idx_lst[1]:ft_idx_lst[2]]
            X_prot_ss = X[:,ft_idx_lst[2]:ft_idx_lst[3]]
            X_pep_2 = X[:,ft_idx_lst[3]:ft_idx_lst[4]]
            X_prot_2 = X[:,ft_idx_lst[4]:ft_idx_lst[5]]
            X_pep_dense = X[:,ft_idx_lst[5]:ft_idx_lst[6]].reshape(X.shape[0],50,3)
            X_prot_dense = X[:,ft_idx_lst[6]:ft_idx_lst[7]].reshape(X.shape[0],589,23)

            pred=model(X_pep,X_prot,X_pep_ss,X_prot_ss,X_pep_2,X_prot_2,X_pep_dense,X_prot_dense)
            preds.extend(pred.detach().cpu().numpy().tolist())
            for row in X:
                seqs.append(info(row))

    assert len(preds) == len(seqs)
    ret = [(seqs[i], preds[i]) for i in range(len(preds))]

    return ret

def pep_seq(row):
    seq = ''
    for i in range(50):
        if row[i].item() == 0:
            break
        seq += rev_map[row[i].item()]
    return seq

def prot_seq(row):
    seq = ''
    for i in range(589):
        if row[50 + i].item() == 0:
            break
        seq += rev_map[row[50 + i].item()]
    return seq

pep_prot_seq = lambda row:f'{pep_seq(row)}, {prot_seq(row)}'

def load_checkpoint(filepath):
    ckpt = torch.load(filepath, map_location='cuda:1')
    model = ckpt['model']
    model.load_state_dict(ckpt['model_state_dict'])
    for parameter in model.parameters():
        parameter.requires_grad=False
    model.eval()
    return model

n_fold = 5

all_preds = {}

with open(f'input.p', 'rb') as f:
    X = pickle.load(f)
load_succ = True

test_X = torch.from_numpy(X)
test_X = test_X.float()
test_X = test_X.to(device)
test_dataset = TensorDataset(test_X)
test_loader = DataLoader(dataset=test_dataset, batch_size=128, shuffle=False)

info = pep_seq

for fold in range(n_fold):
    model_ckpt = load_checkpoint(f'./model_file/model_full_ckpts_{fold}.pkl') # for inference, can refer to these few lines.
    model_ckpt=model_ckpt.to(device)
    print(f'infering using fold{fold+1}')
    preds = infer(model_ckpt, test_loader, device, info)
    for seq, pred in preds:
        if seq not in all_preds.keys():
            all_preds[seq] = []
        all_preds[seq].append(pred)

with open(f'pred.p', '+wb') as f:
    pickle.dump(all_preds, f)